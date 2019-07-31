import unittest
from dimepy import Scan

class ScanTest(unittest.TestCase):
    def test_scan(self):
        # Check if the Scan object fails without pymzml scan
        self.assertRaises(AttributeError, Scan, None, None)

if __name__ == '__main__':
    unittest.main()
