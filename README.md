# Scipy2020
MoSDeF Poster for Scipy2020

## Submitted

### Keyword 
- mosdef, mbuild, foyer, molecular-dynamics, workflow, system-initialization, reproducibility, true

### Abstract
Reproducibility is a prominent issue for many scientific disciplines. 
This issue exists not only in experimental studies but also in computational research. 
Difficulties reproducing data from experiments can be attributed to equipment accuracy or the intrinsic complexity of these studies, i.e. biological systems.
However, in computational research, one would expect calculations to be reproducible with great accuracy, either to machine precision or at least to the precision that parallel processing can provide. 
Yet, many computational studies are difficult to reproduce with sufficient accuracy, if at all.
There are many common reasons, including, but not limited to, software with restrictive licenses that prevent auditing of the source code, a lack of version control and poor documentation, and computational workflows that depend largely on manual input. 
One possible solution to alleviate many of these issues in computational research is to develop tools that automate most/all steps of a computational workflow. In brief, a typical computational chemistry workflow is composed of three steps: system preparation, simulation, and analysis. 
The two latter steps, running and analyzing the simulation, already have mature, open source, and performant libraries with thriving communities built around them. 
System initialization, however, lacks the necessary software infrastructure and is a frequent bottleneck and source of errors for practitioners. 
To this end, we have been developing the Molecular Simulation Design Framework (MoSDeF), a Python-based suite of libraries to help automate and mitigate the fallibility of this step in the computational chemistry workflow. 
By providing researchers with systematic methods and libraries to initialize their simulation systems, we make this step more scriptable, automatable, and hence, reproducible. 
By hosting the open and version-controlled source code on GitHub (https://github.com/mosdef-hub), we encourage users to interact and contribute to the project in a transparent manner. 
Up to this point, much work has been done to build the ground for the project. 
However, there are still many existing obstacles that need to be overcome before MoSDeF can take off as a true open source project. 
The most prominent issue is to attract early adopters that could contribute and shape the future development of this project. 
Due to the scope of this project and the size of the computational chemistry field, the intended user base is relatively small, largely comprised of researchers who may lack software development experience and may be unfamiliar with current best practices. 
This creates a significant barrier for contribution to this project, more than we initially expected. 
Moreover, many computational research groups tend to stick with their “legacy”/ad hoc code, which have been used and passed around their group for years or possibly decades. 
Early efforts to resolve these issues involve grants from the National Science Foundation, with the goal of fostering a community to start using, developing, and stress-testing these tools. 
The overall goal is to shift the field from using closed-source, ad hoc software to open and collaborative software that enables reproducible computational research.

### Guidelines

A virtual conference is a new undertaking for all of us, and we're excited to see what kind of "posters" you come up with!  Without the confines of a physical presentation session, please don't feel obligated to stick to traditional poster formats. While you are welcome to create a single PDF, you are at liberty to create an entire website.  Have a cool new viz tool?  Show us some interactive plots!  New library?  How about a Binder demo showing us what it can do! Consider including a short video providing additional context to your “poster” (2 to 5 minutes).

Your poster - whether it is a single PDF file or interative webpage - must be hosted externally. Please email the link to scipy@enthought.com by June 26th, and we will use this link on the SciPy website.  Feel free to contact us if you have any questions.

The conference is also pleased to announce a best poster contest. We will ask all registered attendees to vote on the poster they deem “best” based on the level of the research, quality of the poster, innovation and clarity of the presentation. The two top winners will receive a certificate and monetary prize.
