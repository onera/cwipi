# Contributing to CWIPI

First of all, thank you for considering contributing to **CWIPI**!


## Development Model

Please note that the core development of **CWIPI** takes place on an internal GitLab server at **ONERA**. The GitHub repository serves as the official mirror and distribution point for the community. All external contributions should follow the workflow described below.


## How to Contribute

### 1. Major Contributions & Feature Proposals
To propose significant features, algorithmic evolutions, or major architectural changes, you *must* open a new GitHub Issue using the **`enhancement`** label to discuss your proposal with the development team. This ensures your work aligns with the project roadmap, and avoids duplicate efforts.

### 2. Development Workflow (Bug Fixes, Features and Documentation)
Whether you are fixing a minor bug, implementing an agreed-upon major feature or improving the documentation, please follow this workflow:

1. Fork the repository on GitHub and clone it locally.
2. Create a new branch branching off from the (up-to-date) **`dev`** branch:
   ```bash
   git checkout dev
   git checkout -b prefix/my-awesome-contribution
    ```
   The prefix (e.g., `bugfix`, `feature` or `docs`) helps to quickly identify the purpose of the branch.
3. Implement your changes locally, ensuring you follow our **Coding Standards** (see below) and write appropriate tests.
4. Push your branch to your GitHub fork.
5. Open a **Pull Request (PR)** from your branch targeting the **`dev`** branch of the official **CWIPI** repository.

**In your Pull Request description, you must clearly specify:**

* The **nature** of the development (what has been changed or added).
* The **motivation** behind it (why this change is necessary or beneficial).
* The related GitHub Issue (if applicable, e.g., `Closes #123`).

### 3. Review Process

Your Pull Request will be reviewed by the core maintainers. Once accepted, it will be integrated into ONERA's internal repository first, tested against our internal CI/CD pipelines, and then pushed back to the public GitHub repository.


## Coding Standards

To maintain a clean and readable codebase, **CWIPI** follows professional standards. While we do not enforce a strict internal manual, our C/C++ code closely adheres to the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).

### Specific Rules:
* **Naming (C):** Use `snake_case` for functions and variables to maintain consistency with standard C libraries.
* **Naming (Python):** Use (upper) `CamelCase` for classes, and `snake_case` for the rest.
* **Indentation:** We use **2 spaces** per indentation level (no tabs).
* **Braces:** The opening brace `{` must be placed **at the end of the line** (not on a new line), as follows:
  ```cpp
  if (condition) {
    // code
  }
  ```

## Attribution
By contributing, you agree that your contributions will be licensed under the project's [LICENSE](LICENSE) (LGPL). All contributors will be acknowledged in the documentation.