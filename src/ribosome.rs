// Copyright 2026 Henrik Mengkrogen
// Developed for the Aarhus University iGEM team.
//
// Licensed under the Apache License, Version 2.0.
// You may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use crate::config::RIBOSOME_SEQUENCE;
use crate::types::RibosomeSimilarity;

pub fn ribosome_similarity(seq_in: &str, n_positions: &[usize]) -> RibosomeSimilarity {
    let seq_bytes = seq_in.as_bytes();
    let ref_bytes = RIBOSOME_SEQUENCE.as_bytes();

    assert!(
        seq_bytes.len() <= ref_bytes.len(),
        "seq_in ({} nt) is longer than the reference ribosome sequence ({} nt)",
        seq_bytes.len(),
        ref_bytes.len()
    );

    let mut mismatched_positions = Vec::new();
    for &pos in n_positions {
        assert!(
            pos < seq_bytes.len(),
            "position {} out of bounds for seq_in",
            pos
        );
        if seq_bytes[pos] != ref_bytes[pos] {
            mismatched_positions.push(pos);
        }
    }

    let total = n_positions.len();
    let matches = total - mismatched_positions.len();
    let identity = if total == 0 {
        1.0
    } else {
        matches as f64 / total as f64
    };

    RibosomeSimilarity {
        matches,
        total,
        identity,
        mismatched_positions,
    }
}
