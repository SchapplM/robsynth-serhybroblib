% Jacobian time derivative of explicit kinematic constraints of
% fourbar2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% WD [5x1]
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)
% [ParkChoPlo1999] Park, FC and Choi, Jihyeon and Ploen, SR: Symbolic formulation of closed chain dynamics in independent coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:22
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function WD = fourbar2DE1_kinconstr_expl_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE1_kinconstr_expl_jacobianD_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar2DE1_kinconstr_expl_jacobianD_mdh_sym_varpar: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE1_kinconstr_expl_jacobianD_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:22:49
% EndTime: 2020-04-24 20:22:49
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [0; 0; 0; 0; 0;];
WD = t1;
