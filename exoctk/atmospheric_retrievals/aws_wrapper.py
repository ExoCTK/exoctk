#!/usr/bin/env python

"""This module serves as a wrapper for executing the ``exoctk``
``atmopsheric_retrievals`` tools using AWS services.

Authors
-------

    - Matthew Bourque

Use
---

    This script is inteneded to be executed via the command line as
    such:

        >>> python aws_wrapper.py

Dependencies
------------

    - boto3
    - paramiko
    - scp
"""

import argparse
import datetime
import getpass
import json
import logging
import os
import socket
import sys
import time

import boto3
import paramiko
from scp import SCPClient


def build_environment(instance):
    """Builds an ``exoctk`` environment on the given AWS EC2 instance

    Parameters
    ----------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    """

    # Establish SSH key
    key = paramiko.RSAKey.from_private_key_file(get_config()['ssh_file'])
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # Connect to the EC2 instance and run commands
    connected = False
    iterations = 0
    while not connected:
        if iterations == 12:
            break
        try:
            client.connect(hostname=instance.public_dns_name, username='ec2-user', pkey=key)
            scp = SCPClient(client.get_transport())
            scp.put('exoctk-aws-build.sh', '~/exoctk-aws-build.sh')
            stdin, stdout, stderr = client.exec_command('chmod 700 exoctk-aws-build.sh && ./exoctk-aws-build.sh')
            connected = True
        except:
            iterations += 1
            time.sleep(5)

    output = stdout.read()

    return output

def configure_logging():
    """Creates a log file that logs the execution of the script

    Returns
    -------
    start_time : obj
        The start time of the script execution
    """

    # Define save location
    log_file = 'logs/aws_wrapper_{}.log'.format(datetime.datetime.now().strftime('%Y-%m-%d-%H-%M'))

    # Make sure no other root lhandlers exist before configuring the logger
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Create the log file
    logging.basicConfig(filename=log_file,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S %p',
                        level=logging.INFO)
    print('Log file initialized to {}'.format(log_file))

    # Log environment information
    logging.info('User: ' + getpass.getuser())
    logging.info('System: ' + socket.gethostname())
    logging.info('Python Version: ' + sys.version.replace('\n', ''))
    logging.info('Python Executable Path: ' + sys.executable)

    start_time = time.time()

    return start_time


def create_ec2():
    """Create an AWS EC2 instance with the given launch template ID
    from the AWS config file.

    Returns
    -------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    """

    ec2 = boto3.resource('ec2')
    LaunchTemplate = {'LaunchTemplateId': get_config()['template_id']}
    instances = ec2.create_instances(
        LaunchTemplate=LaunchTemplate,
        MaxCount=1,
        MinCount=1)
    instance = instances[0]

    logging.info('Launched EC2 instance {}'.format(instance.id))

    instance.wait_until_running()
    instance.load()

    return instance


def get_config():
    """Return a dictionary that holds the contents of the ``jwql``
    config file.

    Returns
    -------
    settings : dict
        A dictionary that holds the contents of the config file.
    """

    config_file_location = 'aws_config.json'

    if not os.path.isfile(config_file_location):
        raise FileNotFoundError('Missing AWS configuration file ("aws_config.json")')

    with open(config_file_location, 'r') as config_file:
        settings = json.load(config_file)

    return settings


def log_execution_time(start_time):
    """Logs the execution time of the script.

    Parameters
    ----------
    start_time : obj
        The start time of the script execution
    """

    end_time = time.time()

    # Log execution time
    hours, remainder_time = divmod(end_time - start_time, 60 * 60)
    minutes, seconds = divmod(remainder_time, 60)
    logging.info('Script Execution Time: {}:{}:{}'.format(int(hours), int(minutes), int(seconds)))


def log_output(output):
    """Logs the standard output of the EC2 instance

    Parameters
    ----------
    output : str
        The standard output of the EC2 instance
    """

    output = output.replace('\t', '  ').replace('\r', '').replace("\'", "").split('\n').decode("utf-8")
    for line in output:
        logging.info(line)


def terminate_ec2(instance):
    """Terminates the given AWS EC2 instance

    Parameters
    ----------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    """

    ec2 = boto3.resource('ec2')
    ec2.instances.filter(InstanceIds=[instance.id]).terminate()

    print('Terminated EC2 instance {}'.format(instance.id))


if __name__ == '__main__':

    start_time = configure_logging()
    instance = create_ec2()
    output = build_environment(instance)
    terminate_ec2(instance)
    log_output(output)
    log_execution_time(start_time)
