% Calculate minimal parameter regressor of potential energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m1DE_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:34
% EndTime: 2020-06-30 17:39:34
% DurationCPUTime: 0.03s
% Computational Cost: add. (23->12), mult. (48->20), div. (0->0), fcn. (48->8), ass. (0->13)
t42 = sin(qJ(1));
t46 = cos(qJ(1));
t38 = g(1) * t46 + g(2) * t42;
t40 = sin(qJ(3));
t44 = cos(qJ(3));
t47 = t40 * g(3) - t38 * t44;
t45 = cos(qJ(2));
t43 = cos(qJ(4));
t41 = sin(qJ(2));
t39 = sin(qJ(4));
t37 = g(1) * t42 - g(2) * t46;
t36 = t44 * g(3) + t40 * t38;
t1 = [0, -t38, t37, 0, 0, 0, 0, 0, g(3) * t41 - t38 * t45, g(3) * t45 + t38 * t41, 0, 0, 0, 0, 0, t41 * t36 + t47 * t45, t36 * t45 - t47 * t41, -t38, 0, t39 * t37 - t43 * t38, t37 * t43 + t39 * t38;];
U_reg = t1;
