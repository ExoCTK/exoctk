#!/usr/bin/env python

"""This module contains various functions for interacting with AWS for
``exoctk`` atmospheric retrievals.

Authors
-------

    - Matthew Bourque

Use
---

    This script is inteneded to be imported and used by other modules,
    for example:

        from aws_tools import create_ec2
        create_ec2()

Dependencies
------------

    - boto3
    - paramiko
    - scp
"""

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


def build_environment(instance, key, client):
    """Builds an ``exoctk`` environment on the given AWS EC2 instance

    Parameters
    ----------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    key : obj
        A ``paramiko.rsakey.RSAKey`` object.
    client : obj
        A ``paramiko.client.SSHClient`` object.
    """

    logging.info('Building ExoCTK environment')

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
    log_output(output)


def configure_logging():
    """Creates a log file that logs the execution of the script

    Returns
    -------
    start_time : obj
        The start time of the script execution
    """

    # Define save location
    log_file = 'logs/aws_wrapper_{}.log'.format(datetime.datetime.now().strftime('%Y-%m-%d-%H-%M'))

    # Create the subdirectory if necessary
    if not os.path.exists('logs/'):
        os.mkdir('logs/')

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
    """Logs the standard output of the EC2 instance.

    Parameters
    ----------
    output : str
        The standard output of the EC2 instance
    """

    output = output.decode("utf-8")
    output = output.replace('\t', '  ').replace('\r', '').replace("\'", "").split('\n')
    for line in output:
        logging.info(line)


def start_ec2(ssh_file, ec2_id):
    """Create a new EC2 instance or start an existing EC2 instance.

    A new EC2 instance will be created if the supplied ``ec2_id`` is an
    EC2 template ID.  An existing EC2 instance will be started if the
    supplied ``ec2_id`` is an ID for an existing EC2 instance.

    Parameters
    ----------
    ssh_file : str
        Relative path to SSH public key to be used by AWS (e.g.
        ``~/.ssh/exoctk.pem``).
    ec2_id : str
        The AWS EC2 template id (e.g. ``lt-021de8b904bc2b728``) or
        instance ID (e.g. ``i-0d0c8ca4ab324b260``).

    Returns
    -------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    key : obj
        A ``paramiko.rsakey.RSAKey`` object.
    client : obj
        A ``paramiko.client.SSHClient`` object.
    """

    ec2 = boto3.resource('ec2')

    # If the given ec2_id is for an EC2 template, then create the EC2 instance
    if ec2_id.split('-')[0] == 'lt':
        LaunchTemplate = {'LaunchTemplateId': ec2_id}
        instances = ec2.create_instances(
            LaunchTemplate=LaunchTemplate,
            MaxCount=1,
            MinCount=1)
        instance = instances[0]
        logging.info('Launched EC2 instance {}'.format(instance.id))

    # If the given ec2_id is for an existing EC2 instance, then start it
    else:
        instance = ec2.Instance(ec2_id)
        instance.start()
        logging.info('Started EC2 instance {}'.format(ec2_id))

    instance.wait_until_running()
    instance.load()

    # Establish SSH key and client
    key = paramiko.RSAKey.from_private_key_file(ssh_file)
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    return instance, key, client


def stop_ec2(ec2_id, instance):
    """Terminates or stops the given AWS EC2 instance.

    The instance is terminated if the supplied ``ec2_id`` is a EC2
    template ID.  The instance is stopped if the supplied ``ec2_id``
    is an ID for a particular EC2 instance.

    Parameters
    ----------
    ec2_id : str
        The AWS EC2 template id (e.g. ``lt-021de8b904bc2b728``) or
        instance ID (e.g. ``i-0d0c8ca4ab324b260``).
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    """

    ec2 = boto3.resource('ec2')

    # If the given ec2_id is for an EC2 template, then terminate the EC2 instance
    if ec2_id.split('-')[0] == 'lt':
        ec2.instances.filter(InstanceIds=[instance.id]).terminate()
        logging.info('Terminated EC2 instance {}'.format(instance.id))

    # If the given ec2_id is for an existing EC2 instance, then stop it
    else:
        instance.stop()
        logging.info('Stopped EC2 instance {}'.format(ec2_id))


def transfer_from_ec2(instance, key, client, filename):
    """Copy files from EC2 user back to the user

    Parameters
    ----------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    key : obj
        A ``paramiko.rsakey.RSAKey`` object.
    client : obj
        A ``paramiko.client.SSHClient`` object.
    filename : str
        The path to the file to transfer
    """

    logging.info('Copying {} from EC2'.format(filename))

    client.connect(hostname=instance.public_dns_name, username='ec2-user', pkey=key)
    scp = SCPClient(client.get_transport())
    scp.get(filename)


def transfer_to_ec2(instance, key, client, filename):
    """Copy parameter file from user to EC2 instance

    Parameters
    ----------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    key : obj
        A ``paramiko.rsakey.RSAKey`` object.
    client : obj
        A ``paramiko.client.SSHClient`` object.
    filename : str
        The path to the file to transfer
    """

    logging.info('Copying {} to EC2'.format(filename))

    connected = False
    iterations = 0
    while not connected:
        if iterations >= 10:
            logging.critical('Could not connect to {}'.format(instance.public_dns_name))
            break
        try:
            client.connect(hostname=instance.public_dns_name, username='ec2-user', pkey=key)
            scp = SCPClient(client.get_transport())
            scp.put(filename)
            connected = True
        except:
            logging.warning('Couldn not connect to {}, retrying.'.format(instance.public_dns_name))
            time.sleep(5)
            iterations += 1

