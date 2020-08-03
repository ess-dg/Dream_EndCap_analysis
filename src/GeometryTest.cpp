#include <gtest/gtest.h>
#include <Geometry.h>

class GeometryTest : public ::testing::Test {
protected:
  Geometry geom;
  CDTReadout CDTRdout;
  CDTGeometry CDTGeom;

  void SetUp() override {}
  void TearDown() override {}
};



// Only valid values are 3 - 6, test a few invalid ones
TEST_F(GeometryTest, InvalidSumo) {
  CDTRdout.anode = 0;
  CDTRdout.cathode = 0;
  CDTRdout.sumo = 2;
  ASSERT_FALSE(geom.calculate(CDTRdout, CDTGeom));
  CDTRdout.sumo = 7;
  ASSERT_FALSE(geom.calculate(CDTRdout, CDTGeom));
}

// Test the valid ones
TEST_F(GeometryTest, ValidSumo) {
  CDTRdout.anode = 0;
  CDTRdout.cathode = 0;

  for (int sumoID = 3; sumoID < 7; sumoID ++) {
    CDTRdout.sumo = sumoID;
    ASSERT_TRUE(geom.calculate(CDTRdout, CDTGeom));
  }
}

TEST_F(GeometryTest, ValidAnodes) {
  const int MaxAnode{63};
  CDTRdout.cathode = 0;
  for (int sumoID = 3; sumoID <=6; sumoID++) {
    CDTRdout.sumo = sumoID;
    for (int anode = 0; anode <= MaxAnode; anode++) {
      CDTRdout.anode = anode;
      ASSERT_TRUE(geom.calculate(CDTRdout, CDTGeom));
    }
    CDTRdout.anode = MaxAnode + 1;
    ASSERT_FALSE(geom.calculate(CDTRdout, CDTGeom));
  }
}

TEST_F(GeometryTest, ValidCathodes) {
  const int MaxCathode{63};
  const int MaxCathodeS6{63};
  CDTRdout.anode = 0;
  for (int sumoID = 3; sumoID <=5; sumoID++) {
    CDTRdout.sumo = sumoID;
    for (int cathode = 0; cathode <= MaxCathode; cathode++) {
      CDTRdout.cathode = cathode;
      ASSERT_TRUE(geom.calculate(CDTRdout, CDTGeom));
    }
    CDTRdout.cathode = MaxCathode + 1;
    ASSERT_FALSE(geom.calculate(CDTRdout, CDTGeom));
  }

  CDTRdout.sumo = 6; // Larger than the rest
  for (int cathode = 0; cathode <= MaxCathodeS6; cathode++) {
    CDTRdout.cathode = cathode;
    ASSERT_TRUE(geom.calculate(CDTRdout, CDTGeom));
  }
  CDTRdout.cathode = MaxCathodeS6 + 1;
  ASSERT_FALSE(geom.calculate(CDTRdout, CDTGeom));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
