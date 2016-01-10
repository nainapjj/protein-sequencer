import unittest
import numpy as np

from .. import motif_analysis as ma
from ...model_generation_p2 import Indexator as ind

class TestMotifClass(unittest.TestCase):

    def test_calculate_volume_from_coords(self):
        """
        Check to ensure that volume calculation function from coords is working as expected

        Used an online calculator at http://www.had2know.com/academics/tetrahedron-volume-4-vertices.html
        to generate accepted values.
        :return:
        """
        coords = [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                  np.array([9.0, 1.0, -3.0])]
        self.assertAlmostEqual(ma.Motif.volumes_from_coords(coords), 34.0)

        coords = [np.array([-1.0, -2.0, -3.0]), np.array([-5.0, -2.0, -3.0]), np.array([-5.0, -10.0, 0.0]),
                  np.array([-9.0, -1.0, 3.0])]
        self.assertAlmostEqual(ma.Motif.volumes_from_coords(coords), 34.0)

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
        self.assertAlmostEqual(motif.calculate_scaled_volume(motif.coords), 13.8804419)

        motif = ma.Motif("ASDF", [np.array([-1.0, -2.0, -3.0]), np.array([-5.0, -2.0, -3.0]),
                                  np.array([-5.0, -10.0, 0.0]), np.array([-9.0, -1.0, 3.0])], [1, 2, 3, 4])
        self.assertAlmostEqual(motif.calculate_scaled_volume(), 13.8804419)
        print("Scaled volume (%f)" % motif.calculate_scaled_volume(motif.coords))

    def test_volume_commutativity(self):
        """
        Sanity check to ensure that the volume functions are commutative
        :return:
        """
        motif_1 = ma.Motif("abcd", [np.array([1.0, 2.0, 3.0]), np.array([5.0, 2.0, 3.0]), np.array([5.0, 10.0, 0.0]),
                               np.array([9.0, 1.0, -3.0])], [1, 2, 4, 8])
        motif_1_identical = ma.Motif("badc", [np.array([5.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]),
                               np.array([9.0, 1.0, -3.0]), np.array([5.0, 10.0, 0.0])], [2, 1, 8, 4])

        print("Volume of motif_1 (%f) vs. motif 1_identical (%f)" % (motif_1.volumes_from_coords(),
                                                                     motif_1_identical.volumes_from_coords()))
        self.assertAlmostEqual(motif_1.volumes_from_coords(), motif_1_identical.volumes_from_coords())

        print("Scaled volume of motif_1 (%f) vs. motif 1_identical (%f)" % (motif_1.calculate_scaled_volume(),
                                                                     motif_1_identical.calculate_scaled_volume()))
        self.assertAlmostEqual(motif_1.calculate_scaled_volume(), motif_1_identical.calculate_scaled_volume())

if __name__ == "__main__":
    unittest.main()

