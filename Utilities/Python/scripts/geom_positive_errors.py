"""
geom_positive_errors.py
This script catches errors produced during FDS setup of bad geometries cases.
The success condition is the existence of the error message. The error message has
its 'ERROR' label replaced by 'SUCCESS' label, thanks to the MISC parameter
POSITIVE_ERROR_TEST=.TRUE. in the bad geometries cases.
The absence of the 'ERROR' label makes current Firebot to consider the run a success.
Original MATLAB script by Emanuele Gissi (6-7-2017)
"""
import os

def check_positive_errors():
    """Check for expected error messages in bad geometry test case error files"""
    fds_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
    err_dir = os.path.join(fds_dir, 'Verification', 'Complex_Geometry', '')
    cases = [
        {
            'file': 'geom_bad_inconsistent_normals.err',
            'error_string': 'Face normals are probably pointing in the wrong direction'
        },
        {
            'file': 'geom_bad_open_surface.err',
            'error_string': 'Open geometry at edge'
        },
        {
            'file': 'geom_bad_non_manifold_edge.err',
            'error_string': 'Non manifold geometry in adjacent faces at edge'
        },
        {
            'file': 'geom_bad_inverted_normals.err',
            'error_string': 'Face normals are probably pointing in the wrong direction'
        },
        {
            'file': 'geom_bad_non_manifold_vert.err',
            'error_string': 'Vertex        6 not connected.'
        }
    ]   
    all_passed = True
    
    # Check existence of error_string in each file
    for n, case in enumerate(cases, start=1):
        infile = os.path.join(err_dir, case['file'])
        errstring = case['error_string']
        
        # Check if file exists
        if not os.path.exists(infile):
            print(f'Error: File {infile} does not exist. Skipping case.')
            all_passed = False
            continue
        
        # Read file and check for error string
        try:
            with open(infile, 'r') as f:
                errfile_content = f.read()
            
            if errstring not in errfile_content:
                print(f'Error: File {infile}:')
                print(f'  Does not contain the following positive error message:')
                print(f'  <{errstring}>')
                all_passed = False
        except Exception as e:
            print(f'Error reading file {infile}: {e}')
            all_passed = False
    return all_passed

if __name__ == '__main__':
    """Main execution block"""
    success = check_positive_errors()

