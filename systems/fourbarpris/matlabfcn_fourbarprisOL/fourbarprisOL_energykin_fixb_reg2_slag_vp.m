% Calculate inertial parameters regressor of fixed base kinetic energy for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbarprisOL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_energykin_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:07
% EndTime: 2020-05-07 09:52:07
% DurationCPUTime: 0.06s
% Computational Cost: add. (1->1), mult. (15->10), div. (0->0), fcn. (0->0), ass. (0->3)
t2 = qJD(1) ^ 2;
t1 = t2 / 0.2e1;
t3 = [0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, qJD(2) * qJD(1), 0, t2 * qJ(2), qJ(2) ^ 2 * t1 + qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, qJD(3) ^ 2 / 0.2e1, 0, 0, 0, 0;];
T_reg = t3;
