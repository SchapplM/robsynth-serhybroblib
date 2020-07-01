% Calculate minimal parameter regressor of fixed base kinetic energy for
% fourbar2DE2
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
% T_reg [1x3]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:59
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar2DE2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE2_energykin_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar2DE2_energykin_fixb_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE2_energykin_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:58:58
% EndTime: 2020-06-26 17:58:59
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0;];
T_reg = t1;
