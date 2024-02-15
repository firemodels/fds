# Contributing

This page provides a starting point for getting involved with the FDS and Smokeview project.  Both novice and experienced developers are also encouraged to review our [Developer Commit Guidelines](https://github.com/firemodels/fds/wiki/Developer-Commit-Guidelines) for details about source code conventions, etc.

## How Can I Best Contribute to FDS?

Historically, FDS development has been the result of the efforts of a small number of people.  With the move towards web based tools for version control and user support, it has become much easier for outside parties to have an opportunity to contribute to FDS.  We, the developers, recognize that this path is a bit of a tightrope.  We never had as much time and resources as we would like for development and user support, so the aid of the user community could be a great benefit.  On the other side is the historical recognition that very few people and organizations are willing or able to work within the existing framework to make a long term concerted effort towards development.  The result is incomplete algorithms which are not validated, too limited in their application, unworkable within the overall FDS framework, or unmaintainable by the core development team. As a result, these efforts never make their way into FDS.  This is of course a frustration both to the developers and to the organization/individual that expended the effort.  The goal of this document is to discuss ways in which interested users can make contributions to the development of FDS that are of lasting value.

## A Few Words on Open Source Development

Open source development can be a great thing.  The concept of a worldwide pool of software developers working in tandem for the greater good is a noble one.  There are hundreds of thousands of open source projects in existence; however, six months from now 99,999 of 100,000 projects will be for all purposes dead in the water.  Why?

 * Too few developers or not enough resources to bring the project to completion
 * Interest in continued development and/or support following a release wanes
 * Disagreements in the direction of the project result in splinter projects being created which dilute the development effort into nothingness

The successful, long lasting open source projects are few (Linux, Mozilla, and OpenOffice are a few).  They have been those with large enough user+developer bases to survive as project members come and go or have committed, long-term institutional support.  

## What can I do to contribute?

### User Support

If the posting logs of the Discussion Forum and the Issue Tracker are examined, you will find that only a few names represent the bulk of responses to user questions and bug fixes.  A very valuable contribution would be for you to spend time on either.  If you know the answer to a question on the forum, post a response.  If you can compile FDS on your own and have a case that gives an error message, compile a debug version and attempt to locate at least the general source of the error before making an issue report.  Time we don't spend on user support is time we have to spend doing development.

### Documentation

This is an area where you can easily make meaningful contributions and in doing so free up time of the core developers.  We attempt to write clear documentation in both the Theory Manual and in the User's Guide.  However, we do not always explain well concepts that we consider "obvious".  As a result documentation may not be as clear or complete as it should be.  If you have struggled with some FDS concept or had to make multiple tries to get an input to work correctly, then help us by writing changes for that portion of the Theory Manual or User's Guide.  A new LaTeX file is preferred, but even pasting text from the pdf file into a word processor and sending changes as, for example, a Word document would help us.  

### Verification and Validation

More work is always needed in the area of V&V.  You can make valuable contributions.  Here's how.  

_Verification_

The goal of the Verification guide is to have fairly simple input files that demonstrate that some particular aspect of FDS has been properly implemented.  For example, that specifying a time dependent injection of mass into a sealed compartment results in the correct mass addition and the correct pressure rise.  Ideally every input should have a simple verification case.  Individually creating a simple case and documenting it takes little time.  To do this for all the FDS inputs becomes very time consuming.  To [Verification_Case_Setup_Example contribute to verification]: 

  * Identify a feature of FDS that currently lacks a verification test case
  * Create a simple test case to test that feature.  These cases should be designed to run quickly (few seconds to a few minutes at most) and they should be setup to create the bare minimum of output needed to determine success or failure
  * Follow the naming and style conventions in the Verification test suite 
  * Document the test case: What feature is being tested, how to evaluate the outputs, and the result of the evaluation.
  * If FDS outputs require post-processing, write a short Matlab script to do that post-processing and to create any plots
 
_Validation_

The goal of the Validation guide is demonstrate the overall predictive performance of FDS, generally by comparing the results of FDS to experimental data or correlations.  To contribute to FDS validation:

  * Identify an aspect of FDS for which there is either no existing validation cases or limited validation cases
  * Identify an experiment that can address the limitation
    * There must be a publicly available test report 
    * The test report must be complete (detailed experimental setup description, assessment of experimental uncertainty, clear identification of measurements made (type and location), and if raw data has been processed a discussion of how that occurred
    * The experimental data (at a minimum that part used for validation) must be able to be included as part of the FDS repository
  * Create FDS input files to simulate the experiment(s). Keep outputs to the minimum needed to perform the validation 
  * Following the approach taking in the Validation Guide, document the experiment and the FDS validation
  * If FDS outputs require post-processing, write a short Matlab script to do that post-processing and create any plots

### Software Development

There are areas of development where you can make valuable contributions.  A good starting place is to read the [FDS Road Map](https://github.com/firemodels/fds/wiki/FDS-Road-Map).  Since FDS is public domain software, you the users are free to do whatever you desire with the FDS source code.  However, if you wish your efforts to become part of a release you will have to abide by the restrictions in the FDS Road Map.  Key among them are avoiding large speed and memory penalties, working within the existing FDS framework, maintaining the broad applicability of FDS, and maintaining the ease of use of FDS, and avoiding the linking to specialized precompiled and/or non-public domain libraries.

Potential developers must keep in mind that successful development is more than just deriving and coding a new algorithm.  The algorithm must be documented, it must undergo verification and validation (the new version of FDS must successfully execute ALL of the V&V suite), it must have a reasonable bang for the buck (i.e. improvement in FDS capability or performance must outweigh any penalties of time and memory usage), and you must be willing to support issue resolution for any bugs related to the new algorithm for a reasonable period of time (six months to a year after release).

What kinds of development is the core development team likely to incorporate into the official release?  

  * New DEVC or CTRL types developed by a user are very likely to incorporated once it is shown the new types do not break existing features.
  * Minor changes to existing inputs that leverage existing code.  A change like allowing a currently constant input parameter to be dependent on a RAMP has a good likelihood of being incorporated once it is shown to successfully run the V&V suite.
  * Changes to the existing algorithms that improve speed / reduce memory usage without changing the fundamental approach.  If no negative impact is seen in executing the V&V suite, and the changes do not make difficult future improvements to FDS, then a good likelihood of being incorporated.
  * New approaches for physical submodels, for example a new method for computing absorptivities for the radiation solver.  If the new model meets the road map restrictions and is a demonstrable improvement in performance and/or accuracy, then it has a good likelihood of being incorporated.
  * *Changes to the fundamental structure of FDS, such as a compressible flow solver, are highly unlikely to be incorporated.*

## How Do I Make a Contribution?

  * Improvements to the wording in guides and answering questions in the discussion forum.  Please do these at any time.
  * Helping to find bugs, volunteer to aid by posting to the Issue Tracker.  
  * New verification cases.  Examine the latest version of the Verification Guide and if your planned contribution is not there, just drop us line to let us know what you plan to do.  
  * New validation cases.  Please contact us in advance.  We will want to verify that your planned contribution is suitable for a validation case.  Not contacting us in advance risks wasted effort if your contribution is deemed unworthy for inclusion.  It should be noted  that we will not automatically refuse a case just because it shows FDS does not perform well.  Though if poor performance is because FDS lacks the capability (compressible flow for example) and there are no plans to ever develop that capability or because the details of the experiment are not known to the level of detailed required to model it we will refuse it.
  * Software development.  For anything other than a new DEVC or CTRL function type we would want to be contacted very early and as well as provided updates during development.  It should be noted that even if we at first indicate the contribution seems worthy of inclusion, that is no guarantee that the end product will be incorporated.
