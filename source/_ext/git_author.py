import os
import subprocess
from sphinx.util import logging

logger = logging.getLogger(__name__)


def _parse_author(email):
    """Return (github_username, None) or (None, email) based on email format."""
    if email.endswith("@users.noreply.github.com"):
        local = email.split("@")[0]
        if "+" in local:
            return local.split("+", 1)[1], None
        return local, None
    return None, email


def _get_first_author(source_dir, docname):
    rst_path = os.path.join(source_dir, docname + ".rst")
    if not os.path.isfile(rst_path):
        return None, None, None

    try:
        result = subprocess.run(
            ["git", "log", "--diff-filter=A", "--format=%an|%ae", "--", rst_path],
            capture_output=True, text=True, timeout=5,
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None, None, None
        last_line = result.stdout.strip().splitlines()[-1]
        name, email = last_line.split("|", 1)
        github_user, plain_email = _parse_author(email.strip())
        return name.strip(), github_user, plain_email
    except Exception as e:
        logger.debug("git_author: failed for %s: %s", docname, e)
        return None, None, None


def add_author_context(app, pagename, templatename, context, doctree):
    name, github_user, email = _get_first_author(app.srcdir, pagename)
    context["first_author_name"] = name
    context["first_author_github"] = github_user
    context["first_author_email"] = email


def setup(app):
    app.connect("html-page-context", add_author_context)
    return {"version": "1.0", "parallel_read_safe": True}
