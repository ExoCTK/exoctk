The ExoCTK team performs a software release before each JWST proposal cycle. We employ the following procedure for creating
a new release:

    1. Create a new branch for changes related to the version release procedure
    2. Update appropriate version numbers in the `setup.py` file,
    3. Update the release notes
    4. Open, review, and merge pull request with the release procedure changes
    5. Create a new tag/release on GitHub
    6. Upload new version of software to PyPI

Detailed instructions for performing a release are given below:

1. Create a new branch for changes related to the version release procedure

Make sure that your local version of the `develop` is up-to-date. A new branch with the naming convention
vx.y.z should be opened off of the `develop`, where vx.y.z is the version number of the release
(e.g. v0.4.1). This branch should be used for the changes described in the rest of this document.

2. Update the version number in `setup.py`

Update all locations of hard-coded version numbers. For example, the VERSION variable in setup.py should be
updated to the new version number, using the x.y.z convention.

3. Update the release notes

In the release notes, write a concise but detailed description of all of the notable changes that have
occurred since the last release. One way to acquire this information is to scroll through the commit history of
the project, and look for commits in which a pull request was merged.

4. Open, review, and merge pull requests with the release procedure changes

Once you've committed the changes from (2), (3), and (4) in your branch, push your branch to GitHub using
the upstream remote, open two pull requests: one that points to the `main` branch, and one
that points to `develop`. Assign reviewers. Either you or the reviewer should eventually merge these pull
requests.

5. Create a new tag/release on GitHub

Once the pull request into the production branch from (5) has been merged, click on the releases button on the
main page of the repository, then hit the "Draft a new release button". The "Tag version" should be the version
number of the release, the "Target" should be the production branch, the "Release title" should (also) be the
version number of the release, and the "Description" should match that of the changelog entry in (4). Once all
of that information is added, hit the big green "Publish" release button.

6. Upload new version of software to PyPI

To upload the new tagged version of the software to PyPI, run the following:

- python setup.py sdist bdist_wheel
- twine upload -u '$<pypi_username>' -p '$<pypi_password>' --repository-url https://upload.pypi.org/legacy/
  --skip-existing dist/*