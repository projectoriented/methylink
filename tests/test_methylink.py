import pytest
from click.testing import CliRunner

from methylink import main

def test_version():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0