# Local test directory
This folder contains a test script intended to do an end-to-end re-training of Metalign (after data and libraries have been setup). This is intended for dev purposes, and to assist with later integration with Travis CI and addressing issue [#22](https://github.com/nlapier2/Metalign/issues/22).

It requires the installation of [bbmap](https://sourceforge.net/projects/bbmap/) to create a mock community for testing.

Please edit the `bbmapRandomReads=...` line in the `retrain_and_test_metalign.sh` to point to the location of your `bbmap/randomreads.sh` file.
