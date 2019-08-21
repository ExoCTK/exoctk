#!/usr/bin/env python

"""This module serves as a test suite for performing atmospheric
retrievals using AWS.  Users must have a ``aws_config.json`` file in
the ``atmospheric_retrievals`` subpackage with appropriate keys.

Executing this script will:

    - Create a AWS EC2 instance from the template defined in the config
      file.
    - Build an ``exoctk`` software environment
    - Run the ``test_atmopsheric_retrievals.py`` module
    - Send result file(s) back to the user
    - Terminate the EC2 instance

Thus, it should be noted that executing this script will result in AWS
costs.

Authors
-------

    - Matthew Bourque

Use
---

    This script is inteneded to be executed via the command line as
    such:

        >>> python test_aws_tools.py

Dependencies
------------

    - scp
"""

import logging

from aws_tools import build_environment
from aws_tools import configure_logging
from aws_tools import create_ec2
from aws_tools import get_config
from aws_tools import log_execution_time
from aws_tools import log_output
from aws_tools import terminate_ec2
from aws_tools import transfer_output_file


def run_tests(instance, key, client):
    """Run atmospheric retrieval unit tests on EC2 instance.

    Parameters
    ----------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    key : obj
        A ``paramiko.rsakey.RSAKey`` object.
    client : obj
        A ``paramiko.client.SSHClient`` object.
    """

    logging.info('Running unit tests')

    command = (
        'export EXOCTK_DATA=""'
        ' && conda activate exoctk-aws'
        ' && cd exoctk/exoctk/tests'
        ' && pytest -s test_atmospheric_retrievals.py'
    )

    # Connect to the EC2 instance and run commands
    client.connect(hostname=instance.public_dns_name, username='ec2-user', pkey=key)
    stdin, stdout, stderr = client.exec_command(command)
    output = stdout.read()
    log_output(output)


if __name__ == '__main__':

    # Get configuration
    ssh_file = get_config()['ssh_file']
    template_id = get_config()['template_id']

    # Configure logging
    start_time = configure_logging()

    # Initialize EC2 instance
    instance, key, client = create_ec2(ssh_file, template_id)

    # Build ExoCTK environment
    build_environment(instance, key, client)

    # Run unit tests
    run_tests(instance, key, client)

    # Transfer output files to user
    transfer_output_file(instance, key, client, 'exoctk/exoctk/tests/BestFit.txt')

    # Terminate EC2 instance
    terminate_ec2(instance)

    # log the execution time
    log_execution_time(start_time)
