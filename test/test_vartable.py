from nose.tools import eq_, ok_
from ngs_doit.types import Bam
# might need the below
# from pyfakefs.fake_filesystem_unittest import Patcher
# from unittest.mock import patch



def func_blah_test1() -> None:
    '''this tests function X for Y '''
    ok_(True, "True is true")
    actual, expected = 3, 3
    eq_(actual, expected, "Actual {0} is not Expected {1}".format(actual, expected))

    # duplicate dirs are impossible anyway, but still using Counter
def func_blah_test2() -> None:
    '''this tests function Z for Y '''
    ok_(True, "True is true")
