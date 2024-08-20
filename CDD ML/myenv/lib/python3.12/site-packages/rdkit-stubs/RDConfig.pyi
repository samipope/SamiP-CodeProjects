"""
 Configuration for the RDKit Python code

"""
from __future__ import annotations
import os as os
import rdkit as rdkit
import sqlite3 as sqlite3
import sys as sys
__all__ = ['ObsoleteCodeError', 'RDCodeDir', 'RDContribDir', 'RDDataDatabase', 'RDDataDir', 'RDDocsDir', 'RDProjDir', 'RDTestDatabase', 'UnimplementedCodeError', 'defaultDBPassword', 'defaultDBUser', 'molViewer', 'os', 'pythonExe', 'pythonTestCommand', 'rdkit', 'rpcTestPort', 'sqlite3', 'sys', 'usePgSQL', 'useSqlLite']
class ObsoleteCodeError(Exception):
    pass
class UnimplementedCodeError(Exception):
    pass
RDCodeDir: str = '/private/var/folders/4_/xyjwmyz1497g3qhhndg003800000gn/T/cirrus-ci-build/build/temp.macosx-11.0-arm64-cpython-312/rdkit_install/lib/python3.12/site-packages/rdkit'
RDContribDir: str = '/private/var/folders/4_/xyjwmyz1497g3qhhndg003800000gn/T/cirrus-ci-build/build/temp.macosx-11.0-arm64-cpython-312/rdkit_install/share/RDKit/Contrib'
RDDataDatabase: str = '/private/var/folders/4_/xyjwmyz1497g3qhhndg003800000gn/T/cirrus-ci-build/build/temp.macosx-11.0-arm64-cpython-312/rdkit_install/share/RDKit/Data/RDData.sqlt'
RDDataDir: str = '/private/var/folders/4_/xyjwmyz1497g3qhhndg003800000gn/T/cirrus-ci-build/build/temp.macosx-11.0-arm64-cpython-312/rdkit_install/share/RDKit/Data'
RDDocsDir: str = '/private/var/folders/4_/xyjwmyz1497g3qhhndg003800000gn/T/cirrus-ci-build/build/temp.macosx-11.0-arm64-cpython-312/rdkit_install/share/RDKit/Docs'
RDProjDir: str = '/private/var/folders/4_/xyjwmyz1497g3qhhndg003800000gn/T/cirrus-ci-build/build/temp.macosx-11.0-arm64-cpython-312/rdkit_install/share/RDKit/Projects'
RDTestDatabase: str = '/private/var/folders/4_/xyjwmyz1497g3qhhndg003800000gn/T/cirrus-ci-build/build/temp.macosx-11.0-arm64-cpython-312/rdkit_install/share/RDKit/Data/RDTests.sqlt'
defaultDBPassword: str = 'masterkey'
defaultDBUser: str = 'sysdba'
molViewer: str = 'PYMOL'
pythonExe: str = '/private/var/folders/4_/xyjwmyz1497g3qhhndg003800000gn/T/cibw-run-cbcjevmk/cp312-macosx_arm64/build/venv/bin/python'
pythonTestCommand: str = 'python'
rpcTestPort: int = 8423
usePgSQL: bool = False
useSqlLite: bool = True
