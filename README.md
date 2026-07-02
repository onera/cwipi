This branch is only here to publish documentation on github.io server.

Follow these steps to do it:

- Copy the generated html directory of each version in `docs/` folder
- Create a file named `.nojekyll` in each version folder. Otherwise, javascript rendering will be deactivated
- Create or update the `docs/index.html` file to make it redirect to the latest doc version, for example:
  ```html
  <head>
    <meta http-equiv='refresh' content='0; URL=1.3.0/index.html'>
  </head>
  ```
- Change hardcoded references in each html file:
  `find . -name '*.html' -exec sed -i 's@href="/coupling/cwipi/@href="/cwipi/@' {} +`
- Add reference to the new version on previous html files:
  `find . -name '*.html' -exec sed -i 's@<dd><a href="/cwipi/1.3.0/">1.3.0</a></dd>@<dd><a href="/cwipi/1.4.0/">1.4.0</a></dd>\n\n          <dd><a href="/cwipi/1.3.0/">1.3.0</a></dd>@' {} +`
