import os
import json

def load_site_settings(setting_file='site_settings/site_settings.json'):
    '''Load the site settings into a dictionary'''
    with open(setting_file, 'r') as f:
        return json.load(f)
