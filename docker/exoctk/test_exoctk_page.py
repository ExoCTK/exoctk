#!/usr/bin/env python
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
        wait = WebDriverWait(driver, 3600)

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
