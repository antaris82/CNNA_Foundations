# CNNA boundary-resolved face-correspondence transport

This package tests whether Double-History gluing induces a nontrivial face correspondence between parent-fan plaquette nets when the correspondence is inferred from boundary/port DtN data rather than from face labels or an imposed Z3 shift.

Primary script:

```bash
python3 test_boundary_resolved_face_correspondence_transport.py --max-level 5 --quick-suite --random-cycles 10 --identical-cycles 10 --outdir boundary_resolved_face_correspondence_out_L5
python3 test_boundary_resolved_face_correspondence_transport.py --max-level 6 --quick-suite --random-cycles 5 --identical-cycles 5 --fingerprints boundary_ports,K_signed --outdir boundary_resolved_face_correspondence_out_L6
```

The primary anti-smuggling fingerprint is `boundary_ports`: it uses DtN boundary/port profiles and response magnitudes, but excludes K/J as matching features. K/J fingerprints are included as audit features only.

No physical complex scalar, no global J, no C*-norm, no GNS construction and no AQFT net are introduced.
