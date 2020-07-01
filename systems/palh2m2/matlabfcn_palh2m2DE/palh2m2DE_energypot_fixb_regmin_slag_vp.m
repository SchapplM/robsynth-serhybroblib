% Calculate minimal parameter regressor of potential energy for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m2DE_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:45
% EndTime: 2020-06-30 17:56:45
% DurationCPUTime: 0.02s
% Computational Cost: add. (18->11), mult. (40->16), div. (0->0), fcn. (40->8), ass. (0->11)
t40 = sin(qJ(1));
t44 = cos(qJ(1));
t36 = g(1) * t44 + g(2) * t40;
t43 = cos(qJ(2));
t42 = cos(qJ(3));
t41 = cos(qJ(4));
t39 = sin(qJ(2));
t38 = sin(qJ(3));
t37 = sin(qJ(4));
t35 = g(1) * t40 - g(2) * t44;
t1 = [0, -t36, t35, 0, 0, 0, 0, 0, -g(3) * t39 - t36 * t43, -t43 * g(3) + t36 * t39, -t36, 0, 0, 0, 0, 0, -g(3) * t38 - t36 * t42, -t42 * g(3) + t36 * t38, -t36, 0, t37 * t35 - t41 * t36, t35 * t41 + t37 * t36;];
U_reg = t1;
