from click.testing import CliRunner

from methylink.main import base


def test_version():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(base, ["--version"])
        assert result.exit_code == 0
