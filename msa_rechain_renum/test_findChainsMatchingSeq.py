import pytest
from findChainsMatchingSeq import *

# def test_main_basic():
#     argv = ['-o', 'test1_out.txt']
#     assert main(argv) == 0

def test_main_1ubq():
    argv = ['1ubq.pdb.gz', '-o', 'test2_out.txt']
    print argv
    assert main(argv) == 0

def test_main_listfile():
    argv = ['-l', 'pdbs.txt', '-o', 'test3_out.txt']
    assert main(argv) == 0

def test_fastas():
    argv = ['-l', 'fastas.txt', '-o', 'test4_out.txt']
    assert main(argv) == 0


# def test_main_config():
#     argv = ['-c', 'config.txt', '-l', 'pdbs.txt']
#     assert main(argv) == 0
#
# def test_main_2hd5():
#     argv = ['/home/hajosep/062014_UBSRD_DB/2HD5.pdb']
#     print argv
#     assert main(argv) == 0
#
# def test_main_listfile_and_pdbs():
#     argv = ['-l', 'pdbs.txt', '-s', 'mysuffix', '1V74.pdb.gz']
#     print argv
#     assert main(argv) == 0


# def test_mytest():
#     with pytest.raises(SystemExit):
#         f()
#
# class TestClass:
#     def test_one(self):
#         x = 'this'
#         assert 'h' in x
#
#     def test_two(self):
#         x = 'hello'
#         assert hasattr(x, 'check')
#
# # tmpdir is a fixture provided by pytest
# # provided fixtures can be listed out with py.test --fixtures
# def test_needsfiles(tmpdir):
#     print tmpdir
#     assert 0
