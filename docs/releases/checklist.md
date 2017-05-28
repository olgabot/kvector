# Release checklist

- [ ] 1. Check that version numbers in `kvector/kvector/__init__.py` and `kvector/setup.py` match
- [ ] 2. Add release notes in `kvector/docs/releases`
- [ ] 3. Copy release notes to `kvector/HISTORY.rst`

- [ ] 4. Convert `README.md` to RST:
```
pandoc --from=markdown_github --to=rst README.md > README.rst
```

- [ ] 4. Check that the `.rst` files are valid using [`rstcheck`](https://pypi.python.org/pypi/rstcheck/0.5) and the built-in Python RST checker:

```
python setup.py check --restructuredtext
rstcheck README.rst
rstcheck HISTORY.rst
```

- [ ] 5. Create an annotated tag for git and push it:

```
git tag -a v0.2.1 -m "v0.2.1 - Release *with* requirements.txt"
git push origin v0.2.1
```

- [ ] 6. (Optional, only if you messed up) and made the tag too early and later made changes that should be added to the tag, remove the remote tag and force re-add the local tag:

```
git push origin :refs/tags/v0.2.1
git tag -fa v0.2.1
git push --tags
```


- [ ] 7. Do a test run of uploading to PyPI test. This assumes you have set up your PyPI configuration according to [this](http://peterdowns.com/posts/first-time-with-pypi.html).
```
python setup.py register -r pypitest
python setup.py sdist upload -r pypitest
```

- [ ] 8. (Optional, if there's an error) in either the `HISTORY.rst` or `README.rst` file, you will get an error. To narrow down this error, look for the corresponding line in `kvector.egg-info`. It won't be exactly the same line because of the header but it will be closer since `PKG-INFO` concatenates your `README.rst` and `HISTORY.rst` files.

```
warning: unexpected indent on line 615 blah blah I'm pypi RST and I'm way too strict
```

- [ ] 9. Check the [PyPI test server](https://testpypi.python.org/pypi) and make sure everything is there and the RST is rendered correctly.
- [ ] 10. Do a test installation in a `conda` environment using the PyPI test server

```
conda create --yes -n kvector_pypi_test_v0.2.1 --file conda_requirements.txt
# Change to a different directory to make sure you're not importing the `kvector` folder
# Use "pushd" instead of "cd" so it's easy to pop back to the kvector folder
pushd $HOME
source activate kvector_pypi_test_v0.2.1
pip install --index-url https://testpypi.python.org/pypi kvector --extra-index-url https://pypi.python.org/simple
```

- [ ] 11. Make sure that the correct `kvector` from the right environment is getting referenced. `which kvector` should have the following output:

```
$ which kvector
/Users/olga/anaconda3/envs/kvector_pypi_v0.2.1/bin/kvector
```

- [ ] 12. Check that the installation was successful. `kvector -h` should have the following output:

```
$ kvector -h
usage: kvector [-h] {index,validate,psi} ...

Calculate "percent-spliced in" (Psi) scores of alternative splicing on a *de
novo*, custom-built splicing index

positional arguments:
  {index,validate,psi}  Sub-commands
    index               Build an index of splicing events using a graph
                        database on your junction reads and an annotation
    validate            Ensure that the splicing events found all have the
                        correct splice sites
    psi                 Calculate "percent spliced-in" (Psi) values using the
                        splicing event index built with "kvector index"

optional arguments:
  -h, --help            show this help message and exit
```

- [ ] 13. Return to the original kvector folder and upload to PyPI:

```
# Return to the kvector directory from the `pushd`
popd
# Now upload to PyPI
python setup.py register -r pypi
python setup.py sdist upload -r pypi
```

Check that it appears correctly on the [PyPI website](https://pypi.python.org/pypi).


- [ ] 14. Upload to Bioconda:

```
cd ..
# If `bioconda-recipes` is not yet cloned, do this step. Otherwise you can
# skip it
git clone https://github.com/bioconda/bioconda-recipes
# Skip to this next line if you've already cloned bioconda-recipes
cd bioconda-recipes/recipes
git checkout -b kvector_v0.2.1
# Remove existing kvector folders
rm -rf kvector
conda skeleton pypi kvector
git add kvector
git commit -m "Updated kvector to v0.2.1"
git push origin kvector_v0.2.1
# Go up two folders to return to the parent folder
cd ../../
```

- [ ] 14. Check that a new pull request exists at [bioconda-recipes](https://github.com/bioconda/bioconda-recipes)

- [ ] 16. Upload to conda-forge

```
# Clone your forked repo of https://github.com/conda-forge/staged-recipes
git clone https://github.com/YeoLab/staged-recipes
cd staged-recipes/recipes
git checkout -b kvector_v0.2.1
conda skeleton pypi kvector
# Get md5 hash from PyPI
git add kvector
git commit -m "Updated kvector to v0.2.1"
```

- [ ] 17. Create a [pull request](https://github.com/conda-forge/staged-recipes/compare/master...YeoLab:master?expand=1) to `conda-forge` using the YeoLab fork.
