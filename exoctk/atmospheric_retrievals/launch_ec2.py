#!/usr/bin/env python

"""
"""

import json
import os
import time

import argparse
import boto3
import paramiko


def create_ec2():
    """Create an AWS EC2 instance with the given launch template ID
    from the AWS config file.
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

    hostname = instance.public_dns_name

    key = paramiko.RSAKey.from_private_key_file(get_config()['ssh_file'])
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    connected = False
    while not connected:
        try:
            client.connect(hostname=hostname, username='ec2-user', pkey=key)
            stdin, stdout, stderr = client.exec_command('ls -l')
            connected = True
        except:
            time.sleep(5)

    print(stdout.read())


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


if __name__ == '__main__':

    create_ec2()
