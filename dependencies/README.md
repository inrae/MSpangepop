This folder contains the workflow environment definition and the Apptainer image recipe.
Changes to `msprime_box.def` trigger the GitLab CI job that builds and pushes an updated
image to the registry. The workflow currently refers to the `latest` image tag.

Lucien
