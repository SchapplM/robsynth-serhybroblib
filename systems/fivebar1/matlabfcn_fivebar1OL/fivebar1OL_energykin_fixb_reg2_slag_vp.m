% Calculate inertial parameters regressor of fixed base kinetic energy for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fivebar1OL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:13
% EndTime: 2020-04-27 06:13:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (8->4), mult. (28->16), div. (0->0), fcn. (4->4), ass. (0->7)
t8 = qJD(3) ^ 2 / 0.2e1;
t7 = qJD(1) ^ 2 / 0.2e1;
t2 = qJD(1) + qJD(2);
t6 = pkin(2) * qJD(1) * t2;
t1 = qJD(3) + qJD(4);
t5 = pkin(3) * qJD(3) * t1;
t3 = [0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2 / 0.2e1, -cos(qJ(2)) * t6, sin(qJ(2)) * t6, 0, pkin(2) ^ 2 * t7, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 / 0.2e1, cos(qJ(4)) * t5, -sin(qJ(4)) * t5, 0, pkin(3) ^ 2 * t8;];
T_reg = t3;
