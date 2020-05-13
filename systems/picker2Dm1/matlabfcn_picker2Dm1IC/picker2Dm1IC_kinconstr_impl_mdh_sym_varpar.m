% Implicit kinematic constraints of
% picker2Dm1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:55
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = picker2Dm1IC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1IC_kinconstr_impl_mdh_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1IC_kinconstr_impl_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:54:42
% EndTime: 2020-05-11 05:54:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (45->33), mult. (24->15), div. (0->0), fcn. (22->20), ass. (0->9)
t5 = qJ(1) + qJ(2);
t4 = qJ(1) + qJ(8);
t3 = qJ(3) + qJ(9);
t2 = pkin(8) + qJ(5);
t1 = qJ(4) + t5;
t8 = pi / 0.2e1;
t7 = cos(qJ(1));
t6 = sin(qJ(1));
t9 = [t7 * pkin(1) + pkin(4) * cos(t1) + pkin(7) + (cos(t5) + sin(qJ(7))) * pkin(3); t6 * pkin(1) + pkin(4) * sin(t1) + (sin(t5) - cos(qJ(7))) * pkin(3); (cos(t4) - cos(pkin(8))) * pkin(5) + (-cos(t2) - t7) * pkin(1); (sin(t4) - sin(pkin(8))) * pkin(5) + (-sin(t2) - t6) * pkin(1); (-cos(t3) - 0.1e1) * pkin(2) + (cos(qJ(3)) - cos(qJ(6))) * pkin(6); -sin(t3) * pkin(2) + (sin(qJ(3)) - sin(qJ(6))) * pkin(6); -qJ(7) + qJ(10) + t8 + t1; -qJ(11) + t2 - t4; qJ(6) + qJ(12) - t3; t8 - qJ(8) - qJ(9) + qJ(2);];
h = t9;
