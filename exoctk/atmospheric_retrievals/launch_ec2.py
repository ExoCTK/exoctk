#!/usr/bin/env python

"""
"""

import json
import os
import time

import argparse
import boto3
import paramiko


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
        if iterations == 30:
            break
        try:
            client.connect(hostname=instance.public_dns_name, username='ec2-user', pkey=key)
            stdin, stdout, stderr = client.exec_command('ls -l')
            connected = True
        except:
            iterations += 1
            time.sleep(5)

    print(stdout.read())


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

    print('Launched EC2 instance {}'.format(instance.id))

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


def terminate_ec2(instance):
    """Terminates the given AWS EC2 instance

    Parameters
    ----------
    instance : obj
        A ``boto3`` AWS EC2 instance object.
    """

    ec2 = boto3.resource('ec2')
    ec2.instances.filter(InstanceIds=[instance.id]).terminate()


if __name__ == '__main__':

    instance = create_ec2()
    build_environment(instance)
    terminate_ec2(instance)
