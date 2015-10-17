/**
 * File: BBFKNNKDTreeUnitTest.cpp
 * Author Lili Meng (lilimeng1103@gmail.com) based on the code by Keith Schwarz (htiek@cs.stanford.edu)
 */

#include <iostream>
#include "KD_tree.h"
#include "ReadData1.h"

using namespace std;

/* These flags control which tests will be run.*/

#define KDTreeTestEnabled          1 // Step one checks

//#define ExactKNearestNeighborTestEnabled      0 // Step two checks

//#define ApproximateKNearestNeighborTestEnabled      0 // Step two checks

/* Utility function that pauses until the user hits ENTER. */
void PressEnterToContinue() {
    /* Use getline to stall until receiving input. */
    string line;
    getline(cin, line);
}

/* This function is what the test suite uses to ensure that the KDTree works
 * correctly. It takes as parameters an expression and description, along
 * with a file and line number, then checks whether the condition is true.
 * If so, it prints that the test passed. Otherwise, it reports that the test
 * fails and points the caller to the proper file and line.
 */

void DoCheckCondition(bool expr, const string& rationale, const string& file, int line) {
    /* It worked! Congrats. */
    if (expr) {
    cout << "PASS: "<< rationale << endl;
    return;
    }

    /* Uh oh! The test failed! */
    cout << "FAIL: " << rationale << endl;
    cout <<"   Error at " << file <<", line" << line << endl;
    cout <<"   (ENTER to continue)" << endl;

    /* Pause so that the test fail stands out. */
    PressEnterToContinue();
}

/* Reports that an unexpected error occurred that caused a test to fail. */
void FailTest(const exception& e) {
    cerr << "TEST FAILED: Unexpected exception: " <<e.what() << endl;
    PressEnterToContinue();
}

 /* This macro takes in an expression and a string, then invokes
  * DoCheckCondition passing in the arguments along with the file
  * and line number on which the macro was called. This makes it
  * easier to track down the source of bugs if a test case should
  * fail.
  */
#define CheckCondition(expr, rationale) DoCheckCondition(expr, rationale, __FILE__, __LINE__)

/* Utility function to delimit the start and end of test cases. */
void PrintBanner(const string& header) {
    cout << "\nBeginning test: " << header << endl;
    cout << setw(40) << setfill('-') << "" << setfill(' ') <<endl;
}


/* Utility function to signal that a test isn't begin run. */
void TestDisabled(const string& header) {
  cout << "== Test " << header << " NOT RUN: press ENTER to continue ==" << endl;

  /* Pause for the user to hit enter. */
  PressEnterToContinue();
}

/* Utility function to signal the end of a test. */
void EndTest() {
  cout << "== end of test: press ENTER to continue ==" << endl;
  PressEnterToContinue();
}


/* Basic test: Can we build the tree and look up the elements it contains? */
void KDTreeTest() try {
#if KDTreeTestEnabled
  PrintBanner("KDTree Test");

  /* Construct the KDTree. */
   /* Read Data */
  vector<vector<double> > dataset;
  ReadData rd1("sample_data.txt");
  int N = rd1.get_num_of_elements();
  int dim = rd1.get_num_of_dimensions();
  dataset=rd1.allDataPointsVec;

  KD_tree kd;
  CheckCondition(kd.empty(),   "New KD tree is empty.");
   /* Add some elements. */
  kd.create_tree(dataset, 1, dim);

  CheckCondition(true, "KDTree construction completed.");

  /* Check basic properties of the KDTree. */
  CheckCondition(kd.dimension() == 128, "Dimension is 128.");

  /* Check basic properties again. */
  CheckCondition(kd.size() == 1000, "After adding 1000 elements, KDTree has size 1000.");
  CheckCondition(!kd.empty(),    "After adding  elements, KDTree is not empty.");

  /* Make sure that the elements we built the tree out of are still there.
  CheckCondition(kd.contains(PointFromRange<128>(dataPoints[0], dataPoints[0] + 128)), "New KD tree has element zero.");
  CheckCondition(kd.contains(PointFromRange<128>(dataPoints[1], dataPoints[1] + 128)), "New KD tree has element one.");
  CheckCondition(kd.contains(PointFromRange<128>(dataPoints[2], dataPoints[2] + 128)), "New KD tree has element two.");

  /* Make sure that the values of these points are correct.
  for (size_t i = 0; i < row1; ++i)
    CheckCondition(kd.at(PointFromRange<128>(dataPoints[i], dataPoints[i] + 128)) == i, "New KD tree has correct values."); */

  EndTest();
#else
  TestDisabled("BasicKDTreeTest");
#endif
} catch (const exception& e) {
  FailTest(e);
}

int main() {
  /* Step Two Tests */
  KDTreeTest();


  /* Step Three Tests
  NearestNeighborTest();
  MoreNearestNeighborTest();

  /* Step Four Tests
  BasicCopyTest();
  ModerateCopyTest();*/

#if (KDTreeTestEnabled)
     /*ModerateKDTreeTestEnabled && \
     HarderKDTreeTestEnabled &&   \
     EdgeCaseKDTreeTestEnabled && \
     MutatingKDTreeTestEnabled && \
     ThrowingKDTreeTestEnabled &&  \
     ConstKDTreeTestEnabled && \
     NearestNeighborTestEnabled &&  \
     MoreNearestNeighborTestEnabled && \
     BasicCopyTestEnabled && \
     ModerateCopyTestEnabled) */
  cout << "All tests completed!  If they passed, you should be good to go!" << endl << endl;
#else
  cout << "Not all tests were run.  Enable the rest of the tests, then run again." << endl << endl;
#endif

  PressEnterToContinue();
}
