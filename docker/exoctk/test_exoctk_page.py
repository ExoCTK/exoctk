#!/usr/bin/env python
"""
Script that uses Selenium to test the exoctk website by navigating to all the sub-pages,
entering values into those pages as needed, hitting "submit", and then making sure that
the results page loads. If the results page is supposed to include any specific download
links, also download those files and check that they exist.

The heart of running a page is the "params" dictionary. This is a dictionary describing
which page to load, what to do there before submitting, how to submit, how to make sure
the results page has loaded, and what (if anything) to do after the results page has 
loaded.

The dictionary has the following form:

{
    'url': "http://localhost:5000/",
        This is a string containing the base URL. Because this script is intended
        primarily to be run on the same docker container that exoctk is running on, it
        always looks at localhost, at the port that docker opens for connections.

    'extension': "contam_visibility",
        This is a string containing the sub-page. The driver will effective load the URL
        URL/EXTENSION

    'pre_steps': [
        {
            "type": "print_head"
        }
    ],
        This is a list of dictionaries of steps to run *before* the form has been 
        submitted. The dictionary format is described in the `run_action()` documentation.


    'submit_id': "calculate_contam_submit",
        The id="VALUE" value of the form submit button. String.

    'submit_done_type': By.CLASS_NAME,
        A value of the "by" submodule. See above for more. Or the string "simple_wait" if
        it is known that the calculation will be finished in a certain time.

    'submit_done_id': "bk-Figure",
        What form element *will* be present after the calculation has finished that is
        *not* present until it has finished. Almost always a string.

    'submit_done_value': value
        The time to wait in seconds if submit_done_type is "simple_wait"

    'post_steps': [
        {
            'type': "download",
            'find_by': By.ID,
            'id': "download_csv",
            'file': 'WASP-18_b_NIS_SUBSTRIP256_visibility.csv'
        }
    ]
        This is a list of dictionaries of steps to run after the form results have loaded.
        The dictionary format is described in the `run_action()` documentation.
"""
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait
from tabulate import tabulate


def print_page_title(url, options, service):
    with webdriver.Firefox(options=options, service=service) as driver:
        driver.get(url)
        print(f"Title: {driver.title}")
        h1_elements = driver.find_elements(By.TAG_NAME, "h1")
        if len(h1_elements) > 0:
            print(f"Top Title: {h1_elements[0].text}")

def print_table(table_array):
    print(f"\n{tabulate(table_array, headers='firstrow', tablefmt='simple')}")

def print_tables(driver, find_by, id):
    tables = driver.find_elements(find_by, id)
    for table in tables:
        table_list = []
        rows = table.find_elements(By.TAG_NAME, "tr")
        for row in rows:
            row_list = []
            cells = row.find_elements(By.TAG_NAME, "td")
            if len(cells) == 0:
                cells = row.find_elements(By.TAG_NAME, "th")
            for cell in cells:
                row_list.append(cell.text)
            table_list.append(row_list)
        print_table(table_list)

def run_action(driver, action, wait_time=1.0, timeout=10.0):
    """
    Runs arbitrary actions (at least if they're defined below)

    Parameters
    ----------
    driver : selenium web driver
    action : dict (see below)
    wait_time : optional float (default 1.0) 
        Time to wait for a download to finish
    timeout : optional float (default 10.0)
        Time to wait for a button to be clickable.

    The action dictionary must have an "action_type" key that can have any of the values
    below:

    - "print_title": No additional arguments. Prints the page title.
    - "print_head": No additional arguments. Prints the text in the first <h1> tag on the
      page.
    - "print_table": Prints all tables with a given ID. Arguments are:
    
        - "find_by": a By value indicating what parameter to search for
        - "id": a string containing the value of the parameter

    - "print_value": Prints the value of a form element. Arguments are:

        - "find_by": a By value
        - "id": The ID to search for

    - "set_text": Set the value of an input element that can be set by `send_keys()`.
      Arguments are:

            - "find_by": a By value
            - "find_value": string containing the value to search for
            - "value": the value to input to the element

    - "resolve_target": Sets the value of the target field and clicks "resolve".
      Arguments are:

        - "target_id": the ID of the form element where the target can be entered
        - "target": the target to resolve
        - "target_check": ID of a form element that will change when the target is resolved

    - "download": Downloads a file and checks for its existence. Arguments are:

        - "find_by": a By value
        - "id": The value to search for
        - "file": The name of the file that will be downloaded
    """
    wait = WebDriverWait(driver, timeout)

    if action["type"] == "print_title":
        # Print the page title
        print(driver.title)
    elif action["type"] == "print_head":
        # Print the text in the first <h1> element in the page
        h1_elements = driver.find_elements(By.TAG_NAME, "h1")
        if len(h1_elements) > 0:
            print(h1_elements[0].text)
    elif action['type'] == "print_table":
        print_tables(driver, action['find_by'], action['id'])
    elif action['type'] == "print_value":
        element = driver.find_element(action['find_by'], action['id'])
        print(f"Element {action['id']} has value {element.get_attribute('value')}")
    elif action["type"] == "set_text":
        # Set a field whose value can be set with `send_keys()`
        element = driver.find_element(action["find_by"], action["find_value"])
        print(f"Set {action['find_value']} from {element.get_attribute('value')} to {action['value']}")
        element.clear()
        element.send_keys(action["value"])
        print(f"Element {action['find_value']} value is {element.get_attribute('value')}")
    elif action["type"] == "resolve_target":
        # Set the target of observation to the provided value, trigger it.
        # Print out the provided check element before and after to make sure it took.
        check_element = driver.find_element(By.ID, action["target_check"])
        check_value = check_element.get_attribute('value')
        print(f"Element {action['target_check']} before target resolution: {check_value}")
        target = driver.find_element(By.ID, action["target_id"])
        target.send_keys(action['target'])
        resolve = wait.until(EC.element_to_be_clickable((By.ID, "resolve_submit")))
        wait.until(EC.invisibility_of_element_located((By.ID, "MathJax_Message")))
        resolve.click()
        check_element = driver.find_element(By.ID, action["target_check"])
        check_value = check_element.get_attribute('value')
        print(f"Element {action['target_check']} after target resolution: {check_value}")
    elif action['type'] == "download":
        element = driver.find_element(action['find_by'], action['id'])
        element.click()
        driver.implicitly_wait(wait_time)
        file_location = Path(f"/root/Downloads/{action['file']}")
        if file_location.is_file():
            print(f"Downloaded {action['file']}")
            file_location.unlink()
        else:
            raise FileNotFoundError(f"Downloaded file {action['file']} not found")

def do_submit(driver, params, calculation_timeout=7200):
    """
    Finds the submit button, waits for it to be pushable, and pushes it.
    """
    submit_element = WebDriverWait(driver, "10").until(
        EC.element_to_be_clickable((By.ID, params["submit_id"]))
    )
    WebDriverWait(driver, "10").until(
        EC.invisibility_of_element_located((By.ID, "MathJax_Message"))
    )
    submit_element.click()
    if params['submit_done_type'] == "simple_wait":
        driver.implicitly_wait(params['submit_done_value'])
    else:
        submit_done = WebDriverWait(driver, calculation_timeout).until(
            EC.presence_of_element_located((params["submit_done_type"],
            params['submit_done_id']))
        )

def do_form(options, service, params):
    """
    Runs all the elements in the pre_steps in order.
    Submits the form.
    Runs all the elements in the post_steps in order.
    """
    run_str = f"Running Test {params['name']}"
    print()
    print(run_str)
    print("-" * len(run_str))
    with webdriver.Firefox(options=options, service=service) as driver:
        # Load groups-integrations page
        driver.get(params['url'] + params['extension'])

        for action in params['pre_steps']:
            run_action(driver, action)

        # Submit the form
        print("Starting Calculation")
        do_submit(driver, params)
        print("Finished Calculation")

        # Do all the post-calculation steps
        for action in params["post_steps"]:
            run_action(driver, action)
    print("-" * len(run_str))
    print()

if __name__ == "__main__":
    # Let's make functions
    options = Options()
    options.binary_location = r'/usr/bin/firefox-esr'
    options.add_argument("--headless")
    service = Service('/usr/bin/geckodriver')

    url = "http://localhost:5000"
    print_page_title(url, options, service)
    print()

    group_integration_params = {
        'name': 'Groups/Integrations',
        'url': "http://localhost:5000/",
        'extension': "groups_integrations",
        'pre_steps': [
            {"type": "print_head"},
            {
                "type": "resolve_target",
                'target_id': "targname",
                'target': "Wasp 18 b",
                'target_check': "kmag",
            },
        ],
        'submit_id': "calculate_submit",
        'submit_done_type': By.ID,
        'submit_done_id': 'myTable',
        'post_steps': [
            {'type': 'print_table', 'find_by': By.ID, 'id': 'myTable'}
        ]
    }
    do_form(options, service, group_integration_params)

    limb_darkening_params = {
        'name': 'Limb Darkening',
        'url': "http://localhost:5000/",
        'extension': "limb_darkening",
        'pre_steps': [
            {"type": "print_head"},
            {
                "type": "resolve_target",
                'target_id': "targname",
                'target': "Wasp 18 b",
                'target_check': "teff",
            },
        ],
        'submit_id': "calculate_submit",
        'submit_done_type': By.CLASS_NAME,
        'submit_done_id': "bk-Figure",
        'post_steps': [
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_coefficients",
                'file': 'ldc_result.csv'
            },
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_spam_coefficients",
                'file': 'spam_result.csv'
            },
        ]
    }
    do_form(options, service, limb_darkening_params)

    fortney_grid_params = {
        'name': 'Fortney Grid',
        'url': "http://localhost:5000/",
        'extension': "fortney",
        'pre_steps': [
            {"type": "print_head"},
            {
                "type": "set_text",
                "find_by": By.ID,
                "find_value": "rstar",
                "value": "2.0"
            },
        ],
        'submit_id': "calculate_submit",
        'submit_done_type': By.CLASS_NAME,
        'submit_done_id': "bk-Figure",
        'post_steps': [
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_file",
                'file': 'fortney.dat'
            },
        ]
    }
    do_form(options, service, fortney_grid_params)

    generic_grid_params = {
        'name': 'Generic Grid',
        'url': "http://localhost:5000/",
        'extension': "generic",
        'pre_steps': [
            {"type": "print_head"},
            {
                "type": "set_text",
                "find_by": By.NAME,
                "find_value": "temperature",
                "value": "700"
            },
            {
                "type": "set_text",
                "find_by": By.NAME,
                "find_value": "gravity",
                "value": "10"
            },
            {
                "type": "set_text",
                "find_by": By.NAME,
                "find_value": "r_star",
                "value": "5.5"
            },
            {
                "type": "set_text",
                "find_by": By.NAME,
                "find_value": "r_planet",
                "value": "2.0"
            },
        ],
        'submit_id': "calculate_submit",
        'submit_done_type': By.CLASS_NAME,
        'submit_done_id': "bk-Figure",
        'post_steps': [
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_file",
                'file': 'generic.dat'
            },
            {'type': 'print_table', 'find_by': By.ID, 'id': 'myTable'},
        ]
    }
    do_form(options, service, generic_grid_params)

    phase_constraint_params = {
        'name': 'Phase Constraint',
        'url': "http://localhost:5000/",
        'extension': "phase_constraint",
        'pre_steps': [
            {"type": "print_head"},
            {
                "type": "resolve_target",
                'target_id': "targname",
                'target': "Wasp 18 b",
                'target_check': "orbital_period",
            },
        ],
        'submit_id': "calculate_submit",
        'submit_done_type': "simple_wait",
        'submit_done_value': 1.0,
        'post_steps': [
            {
                'type': "print_value",
                'find_by': By.ID,
                'id': "minimum_phase",
            },
            {
                'type': "print_value",
                'find_by': By.ID,
                'id': "maximum_phase",
            },
        ]
    }
    do_form(options, service, phase_constraint_params)

    contam_visibility_only_params = {
        'name': 'Visibility Check',
        'url': "http://localhost:5000/",
        'extension': "contam_visibility",
        'pre_steps': [
            {"type": "print_head"},
            {
                "type": "resolve_target",
                'target_id': "targname",
                'target': "Wasp 18 b",
                'target_check': "ra",
            },
        ],
        'submit_id': "calculate_submit",
        'submit_done_type': By.CLASS_NAME,
        'submit_done_id': "bk-Figure",
        'post_steps': [
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_csv",
                'file': 'WASP-18_b_NIS_SUBSTRIP256_visibility.csv'
            }
        ]
    }
    do_form(options, service, contam_visibility_only_params)

    contam_overlap_single_params = {
        'name': 'Contamination/Overlap (single PA check)',
        'url': "http://localhost:5000/",
        'extension': "contam_visibility",
        'pre_steps': [
            {"type": "print_head"},
            {
                "type": "resolve_target",
                'target_id': "targname",
                'target': "Wasp 18 b",
                'target_check': "ra",
            },
            {"type": "set_text", "find_by": By.ID, "find_value": "v3pa", "value": "330"}
        ],
        'submit_id': "calculate_contam_submit",
        'submit_done_type': By.CLASS_NAME,
        'submit_done_id': "bk-Figure",
        'post_steps': [
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_csv",
                'file': 'WASP-18_b_NIS_SUBSTRIP256_visibility.csv'
            }
        ]
    }
    do_form(options, service, contam_overlap_single_params)

    contam_overlap_full_params = {
        'name': 'Contamination/Overlap (full field check)',
        'url': "http://localhost:5000/",
        'extension': "contam_visibility",
        'pre_steps': [
            {"type": "print_head"},
            {
                "type": "resolve_target",
                'target_id': "targname",
                'target': "Wasp 18 b",
                'target_check': "ra",
            },
        ],
        'submit_id': "calculate_contam_submit",
        'submit_done_type': By.CLASS_NAME,
        'submit_done_id': "bk-Figure",
        'post_steps': [
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_csv",
                'file': 'WASP-18_b_NIS_SUBSTRIP256_visibility.csv'
            }
        ]
    }
    do_form(options, service, contam_overlap_full_params)
