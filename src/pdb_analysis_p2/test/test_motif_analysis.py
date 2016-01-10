import unittest
import numpy as np
import math

from .. import motif_analysis as ma
from ...model_generation_p2 import Indexator as ind

class TestMotifClass(unittest.TestCase):
    def test_calculate_area_from_coords(self):
        """
        Check to ensure that area calculation function from coords is working as expected

        Used an online calculator at http://www.had2know.com/academics/triangle-area-perimeter-angle-3-coordinates.html
        to generate accepted values.
        :return:
        """
        coords = [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                  np.array([9.0, 1.0, -3.0])]
        self.assertAlmostEqual(round(ma.Motif.area_from_coords(coords), 3), 17.088)

        coords = [np.array([-1.0, -2.0, -3.0]), np.array([-5.0, -2.0, -3.0]), np.array([-5.0, -10.0, 0.0]),
                  np.array([-9.0, -1.0, 3.0])]
        self.assertAlmostEqual(round(ma.Motif.area_from_coords(coords), 3), 17.088)

        coords = [np.array([1.0, 2.0, 4.0]), np.array([5.0, 2.0, 3.0]), np.array([8.0, 10.0, 0.0]),
                  np.array([9.0, 1.0, -3.0])]
        self.assertAlmostEqual(round(ma.Motif.area_from_coords(coords), 3), 17.095)

    def test_calculate_area_from_lengths(self):
        """
        Check to ensure that Heron's formula is working as expected.

        :return:
        """
        lengths = [3, 4, 5]
        self.assertEqual(ma.Motif.area_from_lengths(lengths), 6.0)

    def test_area_commutativity(self):
        """
        Sanity check to ensure that the volume functions are commutative
        :return:
        """
        motif_1 = ma.Motif("abc", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                                    ], [1, 2, 4])
        motif_1_identical = ma.Motif("cab", [np.array([5.0, 10.0, 0.0]), np.array([1.0, 2.0, 3.0]),
                                              np.array([5.0, 2.0, 3.0])], [4, 1, 2])

        print("Volume of motif_1 (%f) vs. motif 1_identical (%f)" % (motif_1.get_measure(),
                                                                     motif_1_identical.get_measure()))
        self.assertAlmostEqual(motif_1.get_measure(), motif_1_identical.get_measure())

        print("Scaled volume of motif_1 (%f) vs. motif 1_identical (%f)" % (motif_1.get_scaled_measure(),
                                                                     motif_1_identical.get_scaled_measure()))
        self.assertAlmostEqual(motif_1.get_scaled_measure(), motif_1_identical.get_scaled_measure())

    def test_calculate_volume_from_coords(self):
        """
        Check to ensure that volume calculation function from coords is working as expected

        Used an online calculator at http://www.had2know.com/academics/tetrahedron-volume-4-vertices.html
        to generate accepted values.
        :return:
        """
        coords = [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                  np.array([9.0, 1.0, -3.0])]
        self.assertAlmostEqual(ma.Motif.volume_from_coords(coords), 34.0)

        coords = [np.array([-1.0, -2.0, -3.0]), np.array([-5.0, -2.0, -3.0]), np.array([-5.0, -10.0, 0.0]),
                  np.array([-9.0, -1.0, 3.0])]
        self.assertAlmostEqual(ma.Motif.volume_from_coords(coords), 34.0)

    def test_calculate_scaling_volume(self):
        """
        Check to ensure that volume calculation function from edge lengths is working as expected.

        :return:
        """

        # Source: Online calculator at http://www.had2know.com/academics/tetrahedron-volume-6-edges.html
        # to generate accepted value.
        indices = [1, 2, 3, 4]
        self.assertAlmostEqual(ma.Motif.volume_scaling_factor(indices), (1.0 / 6))

        # Tests to see if the new scaling function is the exact same as Stuart's old one.
        indices = [1, 2, 4, 8]
        self.assertAlmostEqual(ma.Motif.volume_scaling_factor(indices), ind.indexator(indices))

        indices = [1, 2, 4, 100]
        self.assertAlmostEqual(ma.Motif.volume_scaling_factor(indices), ind.indexator(indices))


    def test_calculate_scaled_volume(self):
        """
        Check to ensure that scaled volume calculation function are working as expected
        :return:
        """
        motif = ma.Motif("ASDF", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                               np.array([9.0, 1.0, -3.0])], [1, 2, 3, 4])
        print("Scaled volume (%f)" % motif.calculate_scaled_volume())
        self.assertAlmostEqual(motif.calculate_scaled_volume(), 34 / (1.0/6))

        motif = ma.Motif("ASDF", [np.array([-1.0, -2.0, -3.0]), np.array([-5.0, -2.0, -3.0]),
                                  np.array([-5.0, -10.0, 0.0]), np.array([-9.0, -1.0, 3.0])], [1, 2, 3, 4])
        self.assertAlmostEqual(motif.calculate_scaled_volume(), 34 / (1.0/6))
        print("Scaled volume (%f)" % motif.calculate_scaled_volume())

    def test_volume_commutativity(self):
        """
        Sanity check to ensure that the volume functions are commutative
        :return:
        """
        motif_1 = ma.Motif("abcd", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                               np.array([9.0, 1.0, -3.0])], [1, 2, 4, 8])
        motif_1_identical = ma.Motif("badc", [np.array([5.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]),
                               np.array([9.0, 1.0, -3.0]), np.array([5.0, 10.0, 0.0])], [2, 1, 8, 4])

        print("Volume of motif_1 (%f) vs. motif 1_identical (%f)" % (motif_1.get_measure(),
                                                                     motif_1_identical.get_measure()))
        self.assertAlmostEqual(motif_1.get_measure(), motif_1_identical.get_measure())

        print("Scaled volume of motif_1 (%f) vs. motif 1_identical (%f)" % (motif_1.get_scaled_measure(),
                                                                     motif_1_identical.get_scaled_measure()))
        self.assertAlmostEqual(motif_1.get_scaled_measure(), motif_1_identical.get_scaled_measure())

    def test_get_measures(self):
        """
        Test the get_measure and get_measure_scaled functions.
        :return:
        """
        motif_triangle = ma.Motif("abc", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]),
                                          np.array([5.0, 10.0, 0.0])], [1, 2, 3])
        motif_tetrahedral = ma.Motif("abcd", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]),
                                              np.array([5.0, 10.0, 0.0]), np.array([9.0, 1.0, -3.0])], [1, 2, 3, 4])
        # Source: http://www.had2know.com/academics/triangle-area-perimeter-angle-3-coordinates.html
        self.assertAlmostEqual(round(motif_triangle.get_measure(), 3), 17.088)
        # Source: http://keisan.casio.com/exec/system/1223267646
        self.assertAlmostEqual(round(motif_triangle.get_scaled_measure(), 3), 17.088 / .5)

        # Source: http://www.had2know.com/academics/tetrahedron-volume-4-vertices.html
        self.assertAlmostEqual(motif_tetrahedral.get_measure(), 34)
        # Source: http://www.had2know.com/academics/tetrahedron-volume-6-edges.html
        self.assertAlmostEqual(motif_tetrahedral.get_scaled_measure(), 34 * 6)


if __name__ == "__main__":
    unittest.main()

