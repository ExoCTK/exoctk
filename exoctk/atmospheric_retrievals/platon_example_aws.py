#! /usr/bin/env python

"""This module provides an example of how to use the
``platon_wrapper`` module in ExoCTK's ``atmopsheric_retrieval``
subpackage, with focus on utilizing the feature for performing the
processing on AWS.

Authors
-------

    - Matthew Bourque

Use
---

    To run the example, one can execute this script via the command
    line as such:

        >>> python platon_example_aws.py

    To run individual examples, one can import the ``example_aws()``
    function and pass as a parameter which example to run.  Available
    examples include:

        from platon_example_aws import example_aws
        example_aws('emcee')
        example_aws('multinest')

Dependencies
------------

    Users must have a ``aws_config.json`` file present within the
    ``atmospheric_retrievals`` subdirectory.  This file must be of
    a valid JSON format and contain two key/value pairs,
    ``ec2_template_id`` and ``ssh_file``, e.g.:

    {
    "ec2_template_id" : "lt-021de8b904bc2b728",
    "ssh_file" : "~/.ssh/my_ssh_key.pem"
    }

    where the ``ec2_template_id`` contains the ID for an EC2 launch
    template, and ``ssh_file`` points to the SSH public key used for
    logging into an AWS account.
"""

from exoctk.atmospheric_retrievals.aws_tools import get_config
from exoctk.atmospheric_retrievals.platon_wrapper import PlatonWrapper


def example_aws(method):
    """Performs an example run of the emcee and multinest retrievals
    using AWS

    Parameters
    ----------
    method : str
        The method to use to perform the atmopsheric retrieval; can
        either be ``multinest`` or ``emcee``
    """

    ssh_file = get_config()['ssh_file']
    ec2_template_id = get_config()['ec2_template_id']

    # Define the fit parameters
    params = {
        'Rs': 1.19,  # Required
        'Mp': 0.73,  # Required
        'Rp': 1.4,  # Required
        'T': 1200.0,  # Required
        'logZ': 0,  # Optional
        'CO_ratio': 0.53,  # Optional
        'log_cloudtop_P': 4,  # Optional
        'log_scatt_factor': 0,  # Optional
        'scatt_slope': 4,  # Optional
        'error_multiple': 1,  # Optional
        'T_star': 6091}  # Optional

    # Initialize the object, set parameters, and perform retreival
    pw = PlatonWrapper()
    pw.set_parameters(params)
    pw.use_aws(ssh_file, ec2_template_id)
    pw.retrieve(method)


if __name__ == '__main__':

    example_aws('multinest')
    example_aws('emcee')
