#!/bin/sh
# Script for converting a name.ipynb file into name.py and name.rst.include
# Run this script in the examples directory, specifying the example name as an argument
# Example usage: ./generate_example.sh WaterOptimization

# Check for pandoc installation
if ! which pandoc > /dev/null 2>&1; then
    echo "Pandoc is not installed. Install pandoc into your system (NOT python environment) for the conversion to .rst, e.g. via 'sudo apt install pandoc'" >&2
    exit 1
fi

# Check an example directory was provided
example_dir="$1"
if [ -z "${example_dir}" ]; then
  echo "Provide an example_dir directory as an argument" >&2
  exit 1
fi
example_name=$(basename "${example_dir}")
echo "Starting example generation for example '${example_name}'"

# Check the target docs example directory exists, if not create it
docs_examples_dir="../doc/source/examples"

# Workaround for a case where we want the link to remain constant
if [ "${example_name}" = "ASECalculator" ]; then
  target_dir="${docs_examples_dir}/AMSCalculator"
else
  target_dir="${docs_examples_dir}/${example_dir}"
fi

if [ ! -d "${target_dir}" ]; then
  mkdir "${target_dir}"
  echo "Generated documentation directory '${target_dir}' as it did not exist already"
fi

# Find notebook file
nb_files=$(find "${example_dir}" -type f -name "*.ipynb" ! -name "*checkpoint.ipynb")
nb_count=$(echo "${nb_files}" | grep -c .)

# Check the number of .ipynb files found and get file name
if [ "${nb_count}" -eq 0 ]; then
    echo "No .ipynb file found in '${example_dir}'" >&2
    exit 1
elif [ "${nb_count}" -gt 1 ]; then
    echo "Multiple .ipynb files found in '${example_dir}'" >&2
    exit 1
fi

nb_path=$(echo "${nb_files}" | head -n 1)
name=$(basename "${nb_path}" .ipynb)
nb_file="${name}.ipynb"
echo "Using ipynb '${name}' for example generation"

echo "Running black formatter for '${nb_file}'"
$AMSBIN/amspython -m black -t py38  -l 120 "${example_dir}/${nb_file}"

# create the .py file
# remove any get_ipython calls (generated if a cell contains for example_dir !amsmovie)
py_file="${name}.py"
$AMSBIN/amspython -m nbconvert --to python --stdout --no-prompt "${example_dir}/${nb_file}" | sed "1s# python# amspython#; /get_ipython/d" > "${example_dir}/${py_file}"

echo "Generated the python file '${example_dir}/${py_file}'"

echo "Running black formatter for '${py_file}'"
$AMSBIN/amspython -m black -t py38  -l 120 "${example_dir}/${py_file}"

# create the .rst file
# do this via markdown as this gives better control over the pandoc conversion e.g. the width of lines for tables
md_file="${name}.md"
rst_file="${name}.rst"
rst_ipynb_file="${name}.ipynb.rst"
$AMSBIN/amspython -m nbconvert --Exporter.preprocessors="nbconvert_utils.PlamsPreprocessor" --to markdown "${example_dir}/${nb_file}"
pandoc --from markdown --to rst --columns=2000 "${example_dir}/${md_file}" -o "${example_dir}/${rst_file}"

# perform some post-manipulation
# - anonymize the paths 
# - convert python to ipython3
# - remove figure captions
if [ "$(uname)" = "Darwin" ]; then
    sed -i '' -e "
    s#/.*/plams/#/path/plams/#g;
    s#code:: python#code:: ipython3#g;
    /^\.\. figure:: / {n;N;N; d;}
    " "${example_dir}/${nb_file}" "${example_dir}/${rst_file}"
else
    sed -i -e "
    s#/.*/plams/#/path/plams/#g;
    s#code:: python#code:: ipython3#g;
    /^\.\. figure:: / {n;N;N; d;}
    " "${example_dir}/${nb_file}" "${example_dir}/${rst_file}"
fi

# convert subsection headings to sub-sub section headings and add the 'Worked Example' subsection heading
cat <<EOF > "${example_dir}/${rst_ipynb_file}"
Worked Example
--------------

EOF

awk '
  /^-{3,}/ {
    dash_count = length($0)
    tildes = "";              #
    for (i = 1; i <= dash_count; i++) {
      tildes = tildes "~"
    }
    print tildes
    next
  }
  { print }
' "${example_dir}/${rst_file}" >> "${example_dir}/${rst_ipynb_file}"

# move the required files over to the doc directory
cp "${example_dir}/${rst_ipynb_file}" "${target_dir}/"
rm "${example_dir}/${md_file}" "${example_dir}/${rst_file}" "${example_dir}/${rst_ipynb_file}"
echo "Generated the ipynb rst file '${target_dir}/${rst_ipynb_file}'"

# move and generated image files
img_dir="${name}_files"
if [ -d "${example_dir}/${img_dir}" ]; then

  # change a line in any generated .svg files to make them show up correctly in Firefox
  svg_files="${example_dir}/${img_dir}/*.svg"
  for svg_file in ${svg_files}; do
    if [ "$svg_file" != "${svg_files}" ]; then  # avoid error when no svgs are found
      if [ "$(uname)" = "Darwin" ]; then
        sed -i '' 's/xmlns:svg/xmlns/' "${svg_file}"
      else
        sed -i 's/xmlns:svg/xmlns/' "${svg_file}"
      fi
    fi
  done

  if [ -d "${target_dir}/${img_dir}" ]; then rm -rf "${target_dir:?}/${img_dir:?}"; fi
  cp -r "${example_dir}/${img_dir}" "${target_dir}/${img_dir}/"
  rm -r "${example_dir:?}/${img_dir:?}"
  echo "Generated the image directory '${target_dir}/${img_dir}'"
fi

# generate the main template file if it does not exits, and the header and footer files
rst_example_file="${example_name}.rst"
rst_header_file="${example_name}.common_header.rst"
rst_footer_file="${example_name}.common_footer.rst"

if [ ! -e "${target_dir}/${rst_example_file}" ]; then
  cp "${docs_examples_dir}/example.template.rst" "${target_dir}/${rst_example_file}"
  echo "Generated documentation file '${target_dir}/${rst_example_file}' as it did not exist already"
fi
cp "${docs_examples_dir}/example.template.common_header.rst" "${target_dir}/${rst_header_file}"
cp "${docs_examples_dir}/example.template.common_footer.rst" "${target_dir}/${rst_footer_file}"

echo "Generated documentation header file '${target_dir}/${rst_header_file}'"
echo "Generated documentation footer file '${target_dir}/${rst_footer_file}'"

# replace the template placeholders with the example name
if [ "$(uname)" = "Darwin" ]; then
    sed -i '' -e "
    s#<example_name>#${example_name}#g;
    s#<example_file_name>#${name}#g;
    " "${target_dir}/${rst_example_file}" "${target_dir}/${rst_header_file}" "${target_dir}/${rst_footer_file}"
else
    sed -i -e "
    s#<example_name>#${example_name}#g;
    s#<example_file_name>#${name}#g;
    " "${target_dir}/${rst_example_file}" "${target_dir}/${rst_header_file}" "${target_dir}/${rst_footer_file}"
fi





