#ifndef PHIST_FORT_TEST_H
#define PHIST_FORT_TEST_H

/* test for equality and print message on failure */
#define EXPECT_EQ(a,b) call parallel_expect((a)==(b),__LINE__,verbose) 
/* test logicals for equality and print message on failure */
#define EXPECT_EQV(a,b) call parallel_expect((a).eqv.(b),__LINE__,verbose) 
/* test for equality and exit on failure */
#define ASSERT_EQ(a,b) call parallel_assert((a)==(b),__LINE__,verbose) 
/* test logicals for equality and exit on failure */
#define ASSERT_EQV(a,b) call parallel_assert((a).eqv.(b),__LINE__,verbose) 

#endif
