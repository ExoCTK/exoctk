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

    'target_id': "targname",
        The ID to search for to find the field where the target star can be set. String.

    'target': "Wasp 18 b",
        The value to set the target star field *to*. String.

    'target_check': "ra",
        The ID of a field to look at that should change when the target is resolved.
        Currently the value of this field isn't checked, but it likely should be. String.

    'pre_steps': [
        {
            "type": "print_head"
        }
    ],
        This is a list of dictionaries of steps to run *before* the target star has been
        resolved. The dictionary has the form:
            {
                "type": 'TYPE',
            }
        where 'TYPE' is a string that can currently take either the value of "print_title",
        which prints out the title, and "print_head", which prints out the text of the
        first <h1> tag found on the page.
        In the future, this dictionary could also take a "value" field, in case there's
        a desire to set the value of a field before resolving the target.

    'resolved_steps': [
         {"type": "set_value", "find_by": By.ID, "find_value": "v3pa", "value": "330"}
    ],
        This is a list of dictionaries of steps to run *after* the target star has been
        resolved. The dictionary has the form:
            {
                "type": 'TYPE',
                "find_by": 'BY'
                "find_value": 'ID_VALUE',
                "value": 'VALUE'
            }
        where 
            - 'TYPE' is a string that can currently only take the value of "set_value",
              which tells the function to set the value of some form element on the page.
            - 'BY' is a value of the selenium.webdriver.common.by submodule. This tells
              the "find" function what it's using to find the element. By.ID and By.NAME
              are the most common, and look for those elements in the HTML element. It is
              also possible to search by class, by contained text, or by something else.
              Note that if multiple page elements have the same value in whatever you're
              searching for, only the first will be returned.
            - 'ID_VALUE' is the value that is being searched for. For example, to search
              a form for an item with id="ra", 'BY' would be By.ID, and 'ID_VALUE' would
              be "ra". String.
            - 'VALUE': The value to set the item *to*. Whatever class that field is 
              expecting.

    'submit_id': "calculate_contam_submit",
        The id="VALUE" value of the form submit button. String.

    'submit_done_type': By.CLASS_NAME,
        A value of the "by" submodule. See above for more.

    'submit_done_id': "bk-Figure",
        What form element *will* be present after the calculation has finished that is
        *not* present until it has finished. Almost always a string.

    'post_steps': [
        {'type': "print_done"},
        {
            'type': "download",
            'find_by': By.ID,
            'id': "download_csv",
            'file': 'WASP-18_b_NIS_SUBSTRIP256_visibility.csv'
        }
    ]
        This is a list of dictionaries of steps to run after the form results have loaded.
        It has the value
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_csv",
                'file': 'WASP-18_b_NIS_SUBSTRIP256_visibility.csv'
            }
        where
            - 'type' is a string giving the type of action. Accepted values are "table",
              which finds and prints all tables corresponding to a particular tag, 
              "print_done", which prints out a message that the calculation is done,
              and "download", which indicates that Selenium should attempt to download a 
              file from the results page.
            - "find_by" is a by value. See above.
            - "id" is the by value being searched for.
            - "file" is only needed for downloads, and holds the name that the file will
              download as.
}


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

def do_resolve_submit(options, service, params):
    with webdriver.Firefox(options=options, service=service) as driver:
        # Load groups-integrations page
        driver.get(params['url'] + params['extension'])
        wait = WebDriverWait(driver, 7200)

        for action in params['pre_steps']:
            if action["type"] == "print_title":
                print(driver.title)
            elif action["type"] == "print_head":
                h1_elements = driver.find_elements(By.TAG_NAME, "h1")
                if len(h1_elements) > 0:
                    print(h1_elements[0].text)

        # Set the target of observation to the provided value, trigger it.
        # Print out the provided check element before and after to make sure it took.
        check_element = driver.find_element(By.ID, params["target_check"])
        check_value = check_element.get_attribute('value')
        print(f"Element {params['target_check']} before target resolution: {check_value}")
        target = driver.find_element(By.ID, params["target_id"])
        target.send_keys(params['target'])
        resolve = wait.until(EC.element_to_be_clickable((By.ID, "resolve_submit")))
        resolve.click()
        check_element = driver.find_element(By.ID, params["target_check"])
        check_value = check_element.get_attribute('value')
        print(f"Element {params['target_check']} after target resolution: {check_value}")

        for action in params["resolved_steps"]:
            if action["type"] == "set_value":
                element = driver.find_element(action["find_by"], action["find_value"])
                element.send_keys(action["value"])
                print(f"Element {action['find_value']} value is {element.get_attribute('value')}")

        # Submit the form
        submit = wait.until(EC.element_to_be_clickable((By.ID, params["submit_id"])))
        submit.click()
        print("Starting Calculation")
        submit_done = wait.until(EC.presence_of_element_located((params["submit_done_type"], params['submit_done_id'])))

        # Do all the post-calculation steps
        for action in params["post_steps"]:
            if action['type'] == "table":
                print_tables(driver, action['find_by'], action['id'])
            elif action['type'] == "print_done":
                print("Finished Calculation")
            elif action['type'] == "download":
                element = driver.find_element(action['find_by'], action['id'])
                element.click()
                driver.implicitly_wait(2)
                file_location = Path(f"/root/Downloads/{action['file']}")
                if file_location.is_file():
                    print(f"Downloaded {action['file']}")
                else:
                    raise FileNotFoundError(f"Downloaded file {action['file']} not found")
                    

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
        'url': "http://localhost:5000/",
        'extension': "groups_integrations",
        'target_id': "targname",
        'target': "Wasp 18 b",
        'target_check': "kmag",
        'pre_steps': [
            {
                "type": "print_head"
            }
        ],
        'resolved_steps': [],
        'submit_id': "calculate_submit",
        'submit_done_type': By.ID,
        'submit_done_id': 'myTable',
        'post_steps': [
            {'type': 'table', 'find_by': By.ID, 'id': 'myTable'}
        ]
    }
    do_resolve_submit(options, service, group_integration_params)
    print()

    contam_overlap_params = {
        'url': "http://localhost:5000/",
        'extension': "contam_visibility",
        'target_id': "targname",
        'target': "Wasp 18 b",
        'target_check': "ra",
        'pre_steps': [
            {
                "type": "print_head"
            }
        ],
        'resolved_steps': [
#            {"type": "set_value", "find_by": By.ID, "find_value": "v3pa", "value": "330"}
        ],
        'submit_id': "calculate_contam_submit",
        'submit_done_type': By.CLASS_NAME,
        'submit_done_id': "bk-Figure",
        'post_steps': [
            {'type': "print_done"},
            {
                'type': "download",
                'find_by': By.ID,
                'id': "download_csv",
                'file': 'WASP-18_b_NIS_SUBSTRIP256_visibility.csv'
            }
        ]
    }
    do_resolve_submit(options, service, contam_overlap_params)
    print()
