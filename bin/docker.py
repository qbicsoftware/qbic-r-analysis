#!/usr/bin/env python
"""
This little helper script triggers Docker builds for every Dockerfile, if the 
version tag is different from the one that is already on DockerHub. 
"""
import os
import sys
import requests
import subprocess
from requests.auth import HTTPBasicAuth

# The DockerHub API URL, replace 'org' and 'repo' accordingly.
DOCKER_API_URL = "https://hub.docker.com/v2/repositories/{org}/{repo}/tags/"
DOCKER_ORGANISATION = "qbicsoftware"

def main():
    """The main entry point of this script.

    Checks recursivly all Dockerfiles from a given path and checks if their version
    is the most recent one on DockerHub. If not, it triggers Docker builds and pushes
    to DockerHub automatically.
    """
    if not os.environ.get('DOCKER_USER') or not os.environ.get('DOCKER_TOKEN'):
        print("Could not fetch Docker credentials from environment.")
        sys.exit(1)

    # Retrieve all Dockerfiles first
    if len(sys.argv) != 2:
        print("You have to provide a path to the Dockerfiles!")
        sys.exit(1)
    file_path=sys.argv[1]
    file_list = []
    for root, dirs, files in os.walk(file_path):
        for f in files:
             if "Dockerfile" in f:
                 file_list.append(os.path.join(root, f))
    
    if not file_list:
        print("No Dockerfile found, nothing to check")
        return
    
    # Get repos and tags from Dockerfile
    local_repotags = tagsfromdocker(file_list)

    # Check the DockerHub API and trigger an image build + push to DockerHub
    # if the local tag is not in semantic conflict with a remote tag.
    for repo, tag, dockerfile in local_repotags:
        repo = repo.lower()
        if tagondockerhub(repo=repo, tag=tag):
            print("Found tag {tag} for repo {repo} on DockerHub, " 
                "no push required.".format(tag=tag, repo=repo))
            continue
        print("Building Docker image for {}:{}".format(repo, tag))
        buildimage(repo, tag, dockerfile)
        print("Pushing image {}:{} to DockerHub".format(repo, tag))
        pushtodocker(repo, tag)

def pushtodocker(repo, tag):
    """Pushes a docker image to DockerHub.
    """
    name = "{orga}/{repo}:{tag}".format(orga=DOCKER_ORGANISATION,
                                        repo=repo.lower(),
                                        tag=tag)
    pwd = os.environ['DOCKER_TOKEN']
    user = os.environ['DOCKER_USER']
    login_cmd = [
        'docker', 'login', '-p', pwd,
        '-u', user
    ]
    subprocess.check_call(login_cmd)

    push_cmd = [
        'docker', 'push', name
    ]
    subprocess.call(push_cmd)

def buildimage(repo, tag, dockerfile):
    """Builds Docker image with specified name and tag
    """
    name = "{orga}/{repo}:{tag}".format(orga=DOCKER_ORGANISATION,
                                        repo=repo.lower(),
                                        tag=tag)
    full_path = os.path.abspath(os.path.dirname(dockerfile))
    docker_cmd = [
        'docker', 'build', 
        '-t', name, '-f', dockerfile, full_path
    ]
    subprocess.check_call(docker_cmd)

def tagsfromdocker(file_list):
    """Extracts repo names and tags from a Dockerfile list

    Returns a list of tuples (reponame, tag)
    """
    repotags = []

    for dockerfile in file_list:
        name, version = "", ""
        content = parsedockerfile(dockerfile)
        for line in content:
            if "LABEL name=" in line:
                name = line.split("=")[1].strip().replace("\"", "")
            if "LABEL version=" in line:
                version = line.split("=")[1].strip().replace("\"", "")
        if name and version:
            repotags.append((name, version, dockerfile))

    return repotags

def parsedockerfile(file):
    with open(file, 'r') as fh: content = fh.readlines()
    return content

def tagondockerhub(repo="", tag=""):
    """Checks, if a given tag is already given for a repo on DockerHub.

    Returns True if tag is already on DockerHub, else False
    """
    if not repo and tag:
        return False
    resp = requests.get(DOCKER_API_URL.format(org=DOCKER_ORGANISATION, repo=repo))
    if resp.status_code != 200:
        print("The API call failed with {code}, could not check the repo {orga}/{repo}"
            .format(code=resp.status_code, orga=DOCKER_ORGANISATION, repo=repo))
        return False
    results = resp.json()["results"]
    remote_tags = [result["name"] for result in results]
    print(remote_tags)
    return tag in remote_tags
            
if __name__ == '__main__':
    main()

