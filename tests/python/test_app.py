from streamlit.testing.v1 import AppTest
from pathlib import Path


def test_default_app_run_has_no_exception():
    app = AppTest.from_file(Path(__file__).parents[2] / "app.py", default_timeout=30)
    app.run()
    assert not app.exception
    assert app.title[0].value == "SURYA Dynamo v3 · exact C++ port"
