#!/usr/bin/env python
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait
from tabulate import tabulate

# Source - https://stackoverflow.com/a/77099069
# Posted by SoSie
# Retrieved 2026-04-09, License - CC BY-SA 4.0

options = Options()
options.binary_location = r'/usr/bin/firefox-esr'
options.add_argument("--headless")

service = Service('/usr/bin/geckodriver')

with webdriver.Firefox(options=options, service=service) as driver:
    driver.get("http://localhost:5000")
    print(f"Main page title: {driver.title}")

with webdriver.Firefox(options=options, service=service) as driver:
    # Load groups-integrations page
    driver.get("http://localhost:5000/groups_integrations")
    kmag = driver.find_element(By.ID, "kmag")
    print(f"Kmag before resolve: {kmag.get_attribute('value')}")

    # Get the "target name" field, set it to "Wasp 18 b", and resolve it
    target = driver.find_element(By.ID, "targname")
    target.send_keys("Wasp 18 b")
    resolve_target = driver.find_element(By.ID, "resolve_submit")
    resolve_target.click()
    kmag = driver.find_element(By.ID, "kmag")
    print(f"Kmag after resolve: {kmag.get_attribute('value')}")

    # Submit the form
    submit = driver.find_element(By.ID, "calculate_submit")
    submit.click()

    # Wait for the form to be done, then print the output table
    wait = WebDriverWait(driver, 1)
    output_table = wait.until(EC.presence_of_element_located((By.ID, "myTable")))
    tables = driver.find_elements(By.ID, "myTable")
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
        print(tabulate(table_list, headers="firstrow", tablefmt="simple"))
        print()
