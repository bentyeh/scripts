#!/usr/bin/python3

import argparse, getpass, subprocess

parser = argparse.ArgumentParser(
    description='SSH port forwarding',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('loc_port', type=int,
    help='local port from which to forward network requests')
parser.add_argument('--host', type=str, default='localhost',
    help='server hostname (relative to <remote>) receiving forwarded network requests')
parser.add_argument('host_port', type=int,
    help='host port receiving forwarded network requests')
parser.add_argument('remote', type=str,
    help='server hostname receiving forwarded network requests')
parser.add_argument('--bind', type=str, default='localhost',
    help='bind the connection to a specific address')
parser.add_argument('--middle_host', type=str,
    help='intermediary server hostname for 2-hop port forwarding')
parser.add_argument('--middle_port', type=int,
    help='intermediary server port for 2-hop port forwarding')
parser.add_argument('--user', type=str, default=None,
    help='username on <remote> (if using direct forwarding) or <middle_host> (if using 2-hop port forwarding)')

args = parser.parse_args()

if args.bind not in (None, ''):
    args.bind += ':'
else:
    args.bind = ''
if args.user is None:
    args.user = getpass.getuser()

hop_args = (args.middle_host is not None, args.middle_port is not None)
if any(hop_args):
    if not all(hop_args):
        raise ValueError('<middle_host> and <middle_port> must be specified at the same time')
    cmd = (f'ssh -f -L {args.bind}{args.loc_port}:localhost:{args.middle_port} {args.user}@{args.middle_host} '
           f'ssh -N -L {args.middle_port}:{args.host}:{args.host_port} {args.remote}')
else:
    cmd = (f'ssh -Nf -L {args.bind}{args.loc_port}:{args.host}:{args.host_port} {args.user}@{args.remote}')

print(f'Executing command: {cmd}')
subprocess.run(cmd, shell=True)
