Run this to update the dependencies box, temporary image for testing

docker build -t msprime_box ./test_folder/ && apptainer build ./test_folder/msprime_box.sif docker-daemon://msprime_box:latest