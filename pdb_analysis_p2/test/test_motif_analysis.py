import unittest
import numpy as np

from .. import motif_analysis as ma

class TestMotifClass(unittest.TestCase):

    def test_calculate_volume(self):
        """
        Check to ensure that volume calculation function are working as expected
        :return:
        """
        motif = ma.Motif("ASDF", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                               np.array([9.0, 1.0, -3.0])], [])
        self.assertAlmostEqual(motif.calculate_volume(), 34.0)

        motif = ma.Motif("ASDF", [np.array([-1.0, -2.0, -3.0]), np.array([-5.0, -2.0, -3.0]),
                                  np.array([-5.0, -10.0, 0.0]), np.array([-9.0, -1.0, 3.0])], [])
        self.assertAlmostEqual(motif.calculate_volume(), 34.0)

    def test_calculate_scaled_volume(self):
        """
        Check to ensure that scaled volume calculation function are working as expected
        :return:
        """
        motif = ma.Motif("ASDF", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                               np.array([9.0, 1.0, -3.0])], [1, 2, 3, 4])
        print("Scaled volume (%f)" % motif.calculate_scaled_volume())
        self.assertAlmostEqual(motif.calculate_scaled_volume(), 13.8804419)

        motif = ma.Motif("ASDF", [np.array([-1.0, -2.0, -3.0]), np.array([-5.0, -2.0, -3.0]),
                                  np.array([-5.0, -10.0, 0.0]), np.array([-9.0, -1.0, 3.0])], [1, 2, 3, 4])
        self.assertAlmostEqual(motif.calculate_scaled_volume(), 13.8804419)
        print("Scaled volume (%f)" % self.motif.calculate_scaled_volume())

    def test_volume_commutativity(self):
        """
        Sanity check to ensure that the volume functions are commutative
        :return:
        """
        motif_1 = ma.Motif("abcd", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                               np.array([9.0, 1.0, -3.0])], [1, 2, 4, 8])
        motif_1_identical = ma.Motif("badc", [np.array([5.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]),
                               np.array([9.0, 1.0, -3.0]), np.array([5.0, 10.0, 0.0])], [2, 1, 8, 4])

        print("Volume of motif_1 (%f) vs. motif 1_identical (%f)" % (motif_1.calculate_volume(),
                                                                     motif_1_identical.calculate_volume()))
        self.assertAlmostEqual(motif_1.calculate_volume(), motif_1_identical.calculate_volume())

        print("Scaled volume of motif_1 (%f) vs. motif 1_identical (%f)" % (motif_1.calculate_scaled_volume(),
                                                                     motif_1_identical.calculate_scaled_volume()))
        self.assertAlmostEqual(motif_1.calculate_scaled_volume(), motif_1_identical.calculate_scaled_volume())

if __name__ == "__main__":
    unittest.main()

