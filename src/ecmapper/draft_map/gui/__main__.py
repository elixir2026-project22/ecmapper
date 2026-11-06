"""python -m ecmapper.draft_map.gui [MODEL [OUT_DIR] | PREP.pkl BRENDA_map_reactions.json BRENDA_map.svg OUT_DIR] ...

With inputs: a new project of them (OUT_DIR: the folder of the built maps); without: the welcome page.
"""

from __future__ import annotations

import sys

from ..project import Project, options_of
from .session import Session, parse_args


def main(argv=None) -> int:
    args, methods = parse_args(argv)
    from qtpy.QtWidgets import QApplication

    from .escher_svg import web_engine_available
    from .window import MainWindow

    web_engine_available()  # QtWebEngine (SVG export of AI layouts) must be imported before the application is made
    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("ecmapper draft map assignment")
    app.setOrganizationName("ecmapper")
    window = MainWindow(Session(args, methods))
    window.show()
    if args.prep is not None or args.model is not None:  # a new project of these inputs; else the welcome page
        project = Project.new(args.model or args.prep, args.subsystems, args.brenda_json, args.brenda_svg,
                              options_of(args), args.assignment)
        project.meta["build_out"] = str(args.out) if args.out else None
        window.start_project(project)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
