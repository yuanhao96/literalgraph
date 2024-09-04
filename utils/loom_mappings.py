#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import needed libraries
import glob
import json
import os
import pickle
import random
import requests
import time

from tqdm import tqdm  # type: ignore
from typing import Any, Dict, IO, List, Set, Tuple, Union


def gets_json_results_from_api_call(url: str, api_key: str) -> Dict:
    """Function makes API requests and returns results as a json file. API documentation can be found here:
    http://data.bioontology.org/documentation. If a 500 HTTP server-side code is return from "status_code" then the
    algorithm pauses for 90-120 seconds before trying the request again. By default, the process sleeps for 5-20
    seconds after each API call.

    Args:
        url: A string containing a URL to be run against an API.
        api_key: A string containing an API key.

    Return:
        A json-formatted file containing API results.

    Raises:
        An exception is raised if a 500 HTTP server-side code is raised.
    """

    response = requests.get(url, headers={'Authorization': 'apikey token=' + api_key})
    time.sleep(random.randint(10, 20))  # ease rate limiting by sleeping for random intervals

    if response.status_code == 500:
        time.sleep(random.randint(90, 120))  # ease rate limiting by sleeping for random intervals
        response = requests.get(url, headers={'Authorization': 'apikey token=' + api_key})

    return json.loads(response.text)


def writes_data_to_file(file_out: str, results: Set[Tuple[str, Any]]) -> None:
    """Function iterates over set of tuples and writes data to text file locally.

    Args:
        file_out: A filepath to write data to.
        results: A set of tuples, where each tuple represents a mapping between two identifiers.

    Returns:
        None.
    """

    print('Writing results to {location}'.format(location=file_out))

    with open(file_out, 'w') as outfile:
        for res in results:
            outfile.write(res[0] + '\t' + res[1] + '\n')

    outfile.close()

    return None


def processes_api_page_results(content: Dict, source2: str) -> Set:
    """Takes a page of API results and processes them to capture and only those mappings that exist between source1
    and source2. The method returns the results as a set of tuples. Between each batch the process sleeps for 90-120
    seconds to ease the burden on the API.

    Args:
        content: A dictionary of API page results.
        source2: A string naming an ontology. 

    Returns:
        unique_edges: A set of tuples, where each tuple is an edge representing a mapping between source1 and source2.
    """

    unique_edges = set()

    if 'collection' not in content.keys():
        raise KeyError('Something went wrong: {}'.format(content['error']))
    else:
        for result in content['collection']:
            if source2 in result['classes'][1]['links']['ontology']:
                source1_class, source2_class = result['classes'][0]['@id'], result['classes'][1]['@id']
                if '.owl' not in source1_class and '.owl' not in source2_class:
                    unique_edges.add((source1_class, source2_class))

        return unique_edges


def extracts_mapping_data(api_key: str, source1: str, source2: str, file_out: str) -> None:
    """Function uses the BioPortal API to retrieve mappings between two sources. The function batch processes the
    results in chunks of 500, writes the data to a temporary directory and then once all batches have been processed,
    the data is concatenated into a single file.

    Args:
        api_key: A string containing a user BiPortal API key.
        source1: A string naming a source ontology that you want to map from.
        source2: A string naming an ontology that you want to map identifiers from source1 to.
        file_out: A filepath to write data to.

    Returns:
        None.
    """

    print('=' * 50 + '\nRetrieving - {src1} - {src2} Mappings\n'.format(src1=source1, src2=source2) + '=' * 50)

    # get the available resources for mappings to source
    ont_source = 'http://data.bioontology.org/ontologies/{source}/mappings/'.format(source=source1)
    api_results = gets_json_results_from_api_call(ont_source, api_key)
    print('Processing {} Pages of Results'.format(api_results['pageCount']))

    # create temp progress directory to store pages
    temp_progress_storage = '/'.join(file_out.split('/')[:-1]) + '/processed_data'
    os.mkdir(temp_progress_storage)

    # batch process api result pages
    total_pages = list(range(1, int(api_results['pageCount']) + 1))
    n = 500 if len(total_pages) > 5000 else 100
    batches = [total_pages[i:i + n] for i in range(0, len(total_pages), n)]

    for batch in range(0, len(batches)):
        print('\nProcessing batch {} of {}'.format(batch + 1, len(batches) + 1))
        page_results = set()
        for page in tqdm(batches[batch]):
            content = gets_json_results_from_api_call(ont_source + '?page={page}'.format(page=page), api_key)
            pickle.dump(content, open(temp_progress_storage + '/page_{}.pkl'.format(str(page)), 'wb'))  # temp store pg
            page_results |= processes_api_page_results(content, source2)
        writes_data_to_file(file_out + '_{batch_num}'.format(batch_num=batch + 1) + '.txt', page_results)
        time.sleep(random.randint(30, 60))  # ease rate limiting by sleeping for random intervals

    return None


def main() -> None:

    # api_key = input('Please provide your BioPortal API Key: ')
    api_key = '2042a6d7-3ce7-4801-9e4d-d4d6ff132c56'
    # source1 = input('Enter ontology source 1: ').upper()
    # source2 = input('Enter ontology source 2: ').upper()
    source1 = 'DOID'
    source2 = 'HP'

    # create temp directory to store batches
    temp_directory = f'./resources/processed_data/temp_{source1}_{source2}'
    try:
        print('Creating a temporary directory to write API results to: {}'.format(temp_directory))
        os.mkdir(temp_directory)
    except FileExistsError:
        check_input = input('There is already a temp directory, should I overwrite it?: Yes/No')
        if check_input.lower() == 'yes':
            os.remove(temp_directory)
            os.mkdir(temp_directory)
        else:
            new_directory = input('Please provide a name for directory to write data to: ')
            os.mkdir('./resources/processed_data/' + new_directory)

    # run program to map identifiers between source1 and source2
    api_write_location = '/{source1}_{source2}_MAP'.format(source1=source1, source2=source2)
    extracts_mapping_data(api_key, source1, source2, temp_directory + api_write_location)

    # concatenate api data stored in temp directory into single file in the temp directory
    with open(temp_directory + '/{}_{}_MAP.txt'.format(source1.upper(), source2.upper()), 'w') as out:
        for filename in tqdm(glob.glob(temp_directory + '/*.txt')):
            for row in list(filter(None, open(filename, 'r').read().split('\n'))):
                source1_class = '_'.join(row.split('\t')[1].split('/')[-2:])
                source2_class = row.split('\t')[0].split('/')[-1]
                out.write(source1_class + '\t' + source2_class + '\n')
    out.close()

    # delete temp progress storage
    os.remove(temp_directory + '/processed_data')


if __name__ == '__main__':
    main()