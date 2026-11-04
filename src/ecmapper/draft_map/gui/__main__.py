"""python -m ecmapper.draft_map.gui [MODEL [OUT_DIR] | PREP.pkl BRENDA_map_reactions.json BRENDA_map.svg OUT_DIR] ..."""

from __future__ import annotations

import sys

from .session import Session, parse_args


def main(argv=None) -> int:
    args, methods = parse_args(argv)
    from qtpy.QtWidgets import QApplication

    from .window import MainWindow

    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("ecmapper draft map assignment")
    window = MainWindow(Session(args, methods))
    window.show()
    if args.prep is None and args.model is None:
        if not window.new_from_model():  # the inputs from the start dialog
            return 0
    else:
        window.start_loading()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
