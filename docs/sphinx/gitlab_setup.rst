*************
Prerequisites
*************

1. Request access to `ASC RE Gitlab`_ through the `ASC RE Gitlab HPC Accounts`_
   page.
2. Set up your profile on `ASC RE Gitlab`_ with ssh keys. You can follow the
   instructions on the `ASC RE Gitlab User Documentation`_:

   * https://re-git.lanl.gov/help/ssh/index.md#generate-an-ssh-key-pair
   * https://re-git.lanl.gov/help/ssh/index.md#add-an-ssh-key-to-your-gitlab-account

***************************************
Clone cpp\_stub into a local repository
***************************************

1. Navigate to the `upstream repository`_.

2. Copy the ssh URL from the blue Gitlab "Clone" button on the
   `upstream repository`_ web page. The URL should look like the following:

   .. code-block:: bash

      ssh://git@re-git.lanl.gov:10022/aea/stub-repositories/tardigrade.git

3. Navigate to your preferred repository directory on your local computer. In a
   terminal, you can follow the example ``sstelmo`` session below

   .. code-block:: bash

      # Start an ssh session to sstelmo.lanl.gov
      $ ssh -X sstelmo.lanl.gov

      # Note the present working directory (pwd) is your home directory
      $ pwd
      /home/<moniker>

      # OPTIONAL: Create a project space repository directory
      $ mkdir -p /projects/$USER/w13repos
      # Change pwd to repository directory
      $ cd /projects/$USER/w13repos
      # Double check pwd is repository parent directory
      $ pwd
      /projects/<moniker>/w13repos

4. Clone the stub repository using the URL copied in step 2.

   .. code-block:: bash

      # Double check pwd is repository parent directory
      $ pwd
      /projects/<moniker>/w13repos

      # Clone the stub repository
      $ git clone ssh://git@re-git.lanl.gov:10022/aea/stub-repositories/tardigrade.git

5. Rename the local repository directory for your project.

   .. code-block:: bash

      # Double check pwd is repository directory
      $ pwd
      /projects/<moniker>/w13repos

      # Observe the stub repo directory name
      $ ls tardigrade -d
      tardigrade

      # Rename the stub repo directory after your project
      $ mv tardigrade my_project

      # Observe that the stub repo directory no longer exists
      $ ls tardigrade -d
      ls: cannot access 'tardigrade': No such file or directory

      # Observe that your project directory exists
      $ ls my_project -d
      my_project

6. Change to your project's local repository directory

   .. code-block:: bash

      # Double check pwd is repository directory
      $ pwd
      /projects/<moniker>/w13repos

      # Change to your project directory
      $ cd my_project

      # Double check pwd is your project directory
      $ pwd
      /projects/<moniker>/w13repos/my_project

********************************
Create a new upstream repository
********************************

1. Navigate to the W-13 `Material Models`_ Gitlab sub-group.

2. Create a new repository by clicking on the blue "New project" button in the
   upper right corner of the sub-group main page.

   .. note::

      If you do not have the "Developer" or "Maintainer" role assigned to you in
      this sub-group, you will not be able to create a new project directly. You can
      request a role change from the `Material Models`_ sub-group owners. Sub-group
      owners may prefer to create a project for you and make you the owner of that
      project. You can check the `Material Models members`_ list for contact
      information.

3. On the "Create new project" page, follow the link for "Create blank project".

   .. note::

      Gitlab offers a feature to create template projects that may make this
      guide much simpler in the future. Contact the ``tardigrade`` developers and `AEA
      Gitlab group`_ owners to discuss progress on simplified repository setup using
      templates.

3. Enter a name for your project in the "Project name" field. Optionally add a
   "project description" and click the blue "Create project" button.

4. Follow the "Push an existing Git repository" instructions at the bottom of
   the new project webpage.

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ git remote rename origin old-origin
      $ git remote add origin ssh://git@re-git.lanl.gov:10022/aea/material-models/dummy.git
      $ git push -u origin --all
      $ git push -u origin --tags

5. Refresh the Gitlab project webpage and verify that the repository code was pushed
   correctly. You should see a list of source files and this Bitbucket parsed
   ``README.rst`` displayed. You can also review the "main" and "dev" branch from
   the left hand side bar "Repository" > "Branches" menu and the Git tags from the
   "Repository" > "Tags" menu.

6. Remove any issue branches from the ``tardigrade`` project on the "Repository" >
   "Branches" menu. You should keep only the "main" and "dev" branches.

7. If everything looks correct on Gitlab project, you can clean up your local
   repository.

   .. warning::

      WARNING: the ``-D`` option FORCE deletes branches. Triple check the
      command and use with caution. If you're uncertain about this step, contact the
      tardigrade developers for help.

   .. code-block:: bash

      # Remove the tardigrade remote
      $ git remote remove old-origin

      # Ensure that you're on the main branch
      $ git checkout main

      # Remove ALL tardigrade branches except main and dev
      $ git branch | grep -v "main\|dev" | xargs git branch -D

***********************************
Update upstream repository settings
***********************************

Gitlab repositories (a.k.a. 'projects') in the `Material Models`_ Gitlab
sub-group inherit permissions and settings from that sub-group. This includes
inherited minimum roles from the parent `AEA Gitlab group members`. These
default permissions and settings provide access to the AEA group runners on W-13
compute servers and minimize the DevOps work required for new Gitlab projects.
For most developers, these inherited repository settings are appropriate and
only a small number of settings must be updated.

1. Click on the gear icon labeled "Settings" in the lower left sidebar of your
   Gitlab project webpage.

2. Click on the "Repository" menu item that appears in the left sidebar

3. From the "Default branch" > "Expand" page, update the default branch from
   "main" to "dev" and click the blue "Save changes" button.

4. From the "Protected branches" > "Expand" page, protect the "main" and "dev"
   branches according to the needs of your project. The recommended settings are:

   * "allowed to merge"

     * main: Maintainers
     * dev: Developers+Maitainers

   * "allowed to push":

     * main: No one 
     * dev: No one

5. From the "Project Information" > "Members" item at the top of the left side
   bar you can add additional permissions by user and UNIX group.

   .. note::

      Minimum project roles are inherited from `AEA Gitlab group`_ and `Material
      Models`_ sub-group.  Individual projects can elevate roles beyond the minimum,
      but cannot reduce roles.

********************
Enable project CI/CD
********************

The ``tardigrade`` project comes pre-configured to perform continuous integration
(CI) and continuous deployment (CD) on W-13's compute server ``sstelmo`` with
testing performed in and deployment to the `W-13 Python Environments`_.

The CI/CD configuration is found in the ``.gitlab-ci.yml`` file. You can read
more about Gitlab CI/CD configuration in the `ASC RE Gitlab User
Documentation`_: https://re-git.lanl.gov/help/ci/README.md.

No project configuration is required for CI/CD of Merge-Requests to or deployment of the
``dev`` branch. As an alternative to full CI/CD configuration, you may remove the
``git`` operations found in the ``CD.sh`` file, for example found using the
``grep`` command as

.. code-block::

   $ pwd
   /projects/<moniker>/w13repos/my_project

   $ grep git CD.sh
       git config user.name "${GITLAB_USER_NAME}"
       git config user.email "${GITLAB_USER_EMAIL}"
       git remote add oauth2-origin https://gitlab-ci-token:${GITLAB_ACCESS_TOKEN}@re-git.lanl.gov/${CI_PROJECT_PATH}.git
       git tag -a ${production_version} -m "production release ${production_version}" || true
       last_merge_hash=$(git log --pretty=format:"%H" --merges -n 2 | tail -n 1)  # Assume last merge was dev->main. Pick previous
       git tag -a ${developer_version} -m "developer release ${developer_version}" ${last_merge_hash} || true
       git push oauth2-origin --tags

You may also simply remove the ``deploy_build`` job entirely from the
``.gitlab-ci.yml`` file, an example job definition is included below, but the
details may change. The key to identifying the deployment job is the ``stage:
deploy`` attribute and shell commands indicating the CD job definition, e.g.
``script: ./CD.sh``.

.. code-block::
   :linenos:

   deploy_build:
     stage: deploy
     variables:
       GIT_STRATEGY: clone
     script: ./CD.sh
     tags:
       - sstelmo-shell-aea
     only:
       - main
       - dev

The ``pages`` job is a special deploy stage job that builds and deploys
documentation to your project's Gitlab Pages, e.g.
https://aea.re-pages.lanl.gov/stub-repositories/tardigrade. This job should be
retained for building and deploying documentation for your project users.

The ``git`` operations automate micro version bumps during main branch
deployment and are not strictly necessary for CI/CD. The ``deploy_build`` job
performs the CD process and is not required for CI, which is performed by the
``test_build`` job.

The only project configuration required to enable the existing Gitlab CI/CD is
to add a project access token. To add a project access with the naming
convention expected by the CI/CD configuration

1. Click on the gear icon labeled "Settings" in the lower left sidebar of your
   Gitlab project webpage.

2. Click on the "Access Tokens" menu item that appears in the left sidebar

3. Enter the *case-sensitive* name ``GITLAB_ACCESS_TOKEN`` in the "Name" field.

4. Check the ``api`` and ``write_repository`` Scope check boxes. Leave the
   remaining check boxes *unchecked*.

5. Click the blue "Create project access token" button.

6. Copy the text in the "Your new project access token" field.

   .. warning::

      When you navigate away from this page, the access token will *NEVER* be
      visible again. If your copy operation fails or if you overwrite the access token
      in your clipboard, you will need to "revoke" the existing access token from the
      "Active project access tokens" table available on the "Access Tokens" webpage
      and create a new access token from scratch.

      It may be helpful to *TEMPORARILY* copy the access token to an
      intermediate text file for steps 7-10. This access token provides write access
      to your project. *DO NOT SAVE THIS ACCESS TOKEN TO A PLAIN TEXT FILE*.

7. Navigate to the "CI/CD" menu item under "Settings" in the left sidebar.

8. Expand the "Variables" section of the "CI/CD" webpage.

9. Click the blue "Add variable" button.

10. Enter ``GITLAB_ACCESS_TOKEN`` in the "Key" field. This variable name is
    case-sensitive.

11. Paste the access token into the "Value" field.

12. Check both the "Protect Variable" and "Mask Variable" check boxes.

    .. warning::

       Failure to check "Protect Variable" will expose your access token to all
       ASC RE Gitlab runners for all CI/CD pipeline executions on all project
       branches. This may inadvertently expose write access to your project on
       future Gitlab mirrored projects, to users who otherwise have no write access, to
       accidental direct pushes on production branches, or on servers not owned by
       W-13.

    .. warning::

       Failure to check "Mask Variable" will expose your access token in plain
       text in all Gitlab project log files on all servers where the CI/CD is
       performed. It will also expose your access token in plain text on the Gitlab
       CI/CD "Varibles" webpage for all users with project roles of Developer or
       greater access.

13. Click the green "Add variable" button.

14. Click on the "Repository" menu item under the "Settings" item in the left
    sidebar.

15. Expand the "Protected branches" section of the "Repository" webpage.

16. Add the project access token, ``GITLAB_ACCESS_TOKEN``, to the "Allowed to
    push" drop down menu of the "main" and "dev" branches.

*******************
Update project name
*******************

.. note::

   The remaining steps are a truncated version of the `Gitlab Flow`_ workflow.
   Critically, these steps will omit the Gitlab issue creation and Gitlab
   Merge-Request (MR) steps. This step-by-step guide will focus on the Git
   operations performed in the your local repository. The Gitlab MR steps are
   described in greater detail in the `Gitlab Flow`_ documentation.

1. Create a branch for your project name updates using your project's branch
   naming conventions if they exist.

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ git checkout -b feature/project-name-updates
      $ git branch
        dev
      * feature/project-name-updates
        main

2. Search for all instances of ``tardigrade``. The list of occurrences will look
   quite long, but we can search and replace with ``sed`` to avoid manual file
   edits. The session below is an example, the exact output may change but the
   commands should work regardless of project re-organization or evolving features.
   The ellipsis indicates truncated output.

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project

      # Recursive, case-insensitive search and count occurrences
      $ grep -ri tardigrade . --exclude-dir={build,.git} | wc -l
      57

      # Recursive, case-insensitive search and display
      $ grep -ri tardigrade . --exclude-dir={build,.git}
      ...

      # Clean list of files with project name
      $ grep -ri tardigrade . --exclude-dir={build,.git} -l
      ./CMakeLists.txt
      ./docs/api.rst
      ./docs/devops.rst
      ./README.md
      ./set_vars.sh
      ./src/cpp/tardigrade.cpp
      ./src/cpp/tardigrade.h
      ./src/cpp/tests/test_tardigrade.cpp

3. Search and replace from command line

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project

      # Replace lower case occurrences in place
      $ sed -i 's/tardigrade/my_project/g' $(grep -ri tardigrade . --exclude-dir={build,.git} -l)
      $ grep -ri tardigrade . --exclude-dir={build,.git} -l
      ./src/cpp/tardigrade.h

      # Replace upper case occurrences in place
      $ sed -i 's/OVERLAP_COUPLING/MY_PROJECT/g' $(grep -ri tardigrade . --exclude-dir={build,.git} -l)

4. Verify no more occurrences of project name ``tardigrade``

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ grep -ri tardigrade . --exclude-dir={build,.git} | wc -l
      0
      $ grep -ri tardigrade . --exclude-dir={build,.git}
      # no stdout to terminal because no occurrences found
      $ grep -ri tardigrade . --exclude-dir={build,.git} -l
      # no stdout to terminal because no files found

5. Search and replace camelCase project name occurrences, e.g. ``overlapCoupling``.

   .. code-block:: bash

      $ grep -r overlapCoupling . --exclude-dir={build,.git}
      ...
      $ sed -i 's/overlapCoupling/myProject/g' $(grep -r overlapCoupling . --exclude-dir={build,.git} -l)
      $ grep -r overlapCoupling . --exclude-dir={build,.git} -l
      # no stdout to terminal because no files found

6. Find files containing the project in their file name

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ find . -type d \( -name .git -o -name build \) -prune -false -o -name "*tardigrade*"
      ./src/cpp/tardigrade.cpp
      ./src/cpp/tardigrade.h
      ./src/cpp/tests/test_tardigrade.cpp

7. Rename files after current project

   .. note::

      The ``rename`` bash command is common, but not ubiquitous, to UNIX-like
      operating systems. It's reasonably ubiquitous on the most common linux
      distributions. You should find it on ``sstelmo``, but probably won't find it on
      macOS.

   .. code-block:: bash

      $ rename tardigrade myproject $(find . -type d \( -name .git -o -name build \) -prune -false -o -name "*tardigrade*")

8. Commit and push your changes to your "remote" or "fork" repository

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      # Add tracked files and message
      $ git commit -a -m "FEAT: replace tardigrade with my_project throughout repository"
      $ git push origin feature/project-name-updates

You can also perform some cleanup in your documentation directory to remove this
walk-through.

From here, the W-13 workflows would return to the Gitlab webpage and submit a
Merge-Request from the ``feature/project-name-updates`` branch of the upstream
repository to the ``dev`` branch of your "Material Models/my_project"
repository. If the ``.gitlab-ci.yml`` file has been kept, the Merge-Request will
automatically begin running the repository build and test job for continuous
integration (CI). No CI/CD configuration is required for Merge-Requests to or
deployment of the ``dev`` branch.

.. note::

   For Merge-Request and CI/CD of the ``main`` branch, see the previous CI/CD
   configuration section in this setup guide.

For continuing development, W-13 workflows recommend that you should keep the
upstream repository production branches, ``dev`` and ``main``, clean from
development work and *NEVER* develop directly on the ``dev`` and ``main``
branches of your local repository. Limit development work to ``feature/thing``
type branches on your local repo and frequently commit changes and push from the
local feature branch back to the upstream repository.

Happy hacking!
