from pathlib import Path

from scm.plams import Molecule

PATH = Path('.') / 'xyz'

m1 = Molecule(PATH / 'RS1.xyz')
m2 = Molecule(PATH / 'RS2.xyz')


def testYES():
    for i in range(3): assert m1.label(i) == m2.label(i)


def testNO():
    for i in range(3,5): assert m1.label(i) != m2.label(i)
