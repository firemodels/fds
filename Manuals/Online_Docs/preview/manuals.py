"""Identity and presentation metadata; scientific content stays in LaTeX."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Manual:
    key: str
    directory: str
    output: str
    title: str
    subtitle: str
    title_es: str
    subtitle_es: str
    entry: str = 'chapter-1.html'

    @property
    def main(self):
        return self.directory + '.tex'

    @property
    def pdf(self):
        return self.directory + '.pdf'


MANUALS = {
    'tech': Manual('tech', 'FDS_Technical_Reference_Guide', 'tech-guide',
                   'Technical Reference Guide', 'Volume 1 · Mathematical Model',
                   'Guía de referencia técnica', 'Volumen 1 · Modelo matemático', 'overview.html'),
    'user': Manual('user', 'FDS_User_Guide', 'user-guide',
                   'User Guide', 'Fire Dynamics Simulator',
                   'Guía del usuario', 'Simulador de dinámica de incendios'),
    'verification': Manual('verification', 'FDS_Verification_Guide', 'verification-guide',
                           'Verification Guide', 'Volume 2 · Verification',
                           'Guía de verificación', 'Volumen 2 · Verificación'),
    'validation': Manual('validation', 'FDS_Validation_Guide', 'validation-guide',
                         'Validation Guide', 'Volume 3 · Validation',
                         'Guía de validación', 'Volumen 3 · Validación', 'overview.html'),
    'config': Manual('config', 'FDS_Config_Management_Plan', 'config-guide',
                     'Configuration Management Plan', 'Volume 4 · Software Quality Assurance',
                     'Plan de gestión de la configuración', 'Volumen 4 · Garantía de calidad del software'),
}

TECH = MANUALS['tech']
