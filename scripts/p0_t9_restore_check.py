"""Read-only health verification of the unchanged restored Nemotron container."""
from pathlib import Path
import json
import subprocess
import urllib.request

OUT = Path(__file__).resolve().parents[1] / 'outputs/p0_t9/nemotron_recovery'

def main():
    before = json.loads((OUT / 'before.json').read_text())
    current = json.loads(subprocess.check_output(['docker', 'inspect', 'nemotron-nano-agent']))[0]
    assert current['State']['Running']
    for key in ['Id', 'Name', 'Image', 'Path', 'Args', 'Mounts']:
        assert current[key] == before[key], key
    for key, value in before['config'].items():
        assert current['Config'][key] == value, key
    for key, value in before['host_config'].items():
        assert current['HostConfig'][key] == value, key
    base = 'http://127.0.0.1:8002'
    health = urllib.request.urlopen(base + '/health', timeout=10).status
    models = json.load(urllib.request.urlopen(base + '/v1/models', timeout=10))
    ids = [entry['id'] for entry in models['data']]
    previous_ids = [entry['id'] for entry in before['models']['data']]
    assert ids == previous_ids
    payload = {'model': ids[0], 'messages': [{'role': 'user', 'content': 'Reply with the single word OK.'}], 'temperature': 0, 'max_tokens': 64, 'chat_template_kwargs': {'enable_thinking': False}}
    request = urllib.request.Request(base + '/v1/chat/completions', data=json.dumps(payload).encode(), headers={'Content-Type': 'application/json'})
    response = json.load(urllib.request.urlopen(request, timeout=120))
    assert response['choices'] and response['choices'][0]['message'].get('content')
    result = {'container_id': current['Id'], 'running': True, 'original_configuration_unchanged': True, 'health_http_status': health, 'model_ids': ids, 'model_ids_match_before_stop': True, 'requested_model_id': 'nemotron_3_nano_omni', 'requested_model_id_present': 'nemotron_3_nano_omni' in ids, 'short_inference_success': True, 'response': response, 'operational_health': 'PASS', 'requested_model_name_check': 'PASS' if 'nemotron_3_nano_omni' in ids else 'FAIL', 'note': 'The pre-existing served name was preserved; no configuration changes were made.'}
    (OUT / 'restoration.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
