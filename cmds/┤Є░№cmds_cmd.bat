@echo off
rd /s /q dist
rd /s /q build
pyinstaller -D cmds_cmd.py
copy data\bond_table.npy dist\cmds_cmd
copy src\config.ini dist\cmds_cmd
pause