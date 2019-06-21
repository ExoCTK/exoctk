#!/usr/bin/env python

"""
"""

import json
import os

import argparse
import boto3


def create_ec2():
    """Create an AWS EC2 instance with the given launch template ID
    from the AWS config file.
    """

    session = boto3.Session(region_name='us-east-1')
    ec2 = session.resource('ec2')
    LaunchTemplate = {'LaunchTemplateId': get_config()['template_id']}
    instance = ec2.create_instances(
        LaunchTemplate=LaunchTemplate,
        MaxCount=1,
        MinCount=1)

    print('Launched EC2 instance {}'.format(instance))


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
