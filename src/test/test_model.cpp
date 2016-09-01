#include <gtest/gtest.h>
#include <glog/logging.h>

TEST(TestModel,Test0) {
	ASSERT_EQ(0,0);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
