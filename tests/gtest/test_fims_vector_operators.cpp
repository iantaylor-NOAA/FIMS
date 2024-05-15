#include "gtest/gtest.h"
#include "common/fims_vector.hpp"

namespace
{

    // Test exp using multiple input values and types
    // Not worth to write many tests when testing thin database wrappers,
    // third-party libraries, or basic variable assignments.

    TEST(Vector, VectorizedOperations)
    {
        fims::Vector<double> a = {1,2,3,4};
        fims::Vector<double> b = {5,6,7,8};
        fims::Vector<double> c = a + b;//c = {6  8 10 12}
        
        
        // Test operator "+"
        EXPECT_EQ(c[0], 6.0);
        EXPECT_EQ(c[1], 8.0);
        EXPECT_EQ(c[2], 10.0);
        EXPECT_EQ(c[3], 12.0);
        
        fims::Vector<double> d = a - b; //d = {-4 -4 -4 -4}
        // Test operator "-"
        EXPECT_EQ(d[0], -4.0);
        EXPECT_EQ(d[1], -4.0);
        EXPECT_EQ(d[2], -4.0);
        EXPECT_EQ(d[3], -4.0);
        
        fims::Vector<double> e = a * b; //d = {5 12 21 32}
        // Test operator "*"
        EXPECT_EQ(e[0], 5.0);
        EXPECT_EQ(e[1], 12.0);
        EXPECT_EQ(e[2], 21.0);
        EXPECT_EQ(e[3], 32.0);
        
        fims::Vector<double> f = a / b; //d = {0.2 0.3333333 0.4285714 0.5}
        // Test operator "/"
        EXPECT_NEAR(f[0], 0.2, 0.0001);
        EXPECT_NEAR(f[1], 0.3333333, 0.0001);
        EXPECT_NEAR(f[2], 0.4285714, 0.0001);
        EXPECT_NEAR(f[3], 0.5, 0.0001);
        
        fims::Vector<double> pi;
        pi = 3.14159265359;
        EXPECT_NEAR(pi, 3.14159265359, 0.000001);
       
    }

}
